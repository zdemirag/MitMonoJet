#include "MitPhysics/Init/interface/ModNames.h"
#include "MitMonoJet/TreeFiller/interface/FillerXlJets.h"
#include "MitAna/DataTree/interface/PFJetCol.h"

#include "MitMonoJet/DataTree/interface/XlSubJet.h"
#include "MitMonoJet/DataTree/interface/XlFatJet.h"
#include "MitCommon/DataFormats/interface/Vect4M.h"
#include "MitCommon/DataFormats/interface/Vect3.h"
#include "MitCommon/DataFormats/interface/Types.h"
#include "MitCommon/MathTools/interface/MathUtils.h"

#include "QjetsPlugin.h"
#include "Qjets.h"

using namespace mithep;

ClassImp(mithep::FillerXlJets)

//--------------------------------------------------------------------------------------------------
FillerXlJets::FillerXlJets(const char *name, const char *title) :
  BaseMod (name,title),
  fIsData (kTRUE),
  fFillVSubJets (kTRUE),
  fFillTopSubJets (kFALSE),
  fNSubDeclustering (kFALSE),
  fBTaggingActive (kFALSE),
  fQGTaggingActive (kTRUE),
  fQGTaggerCHS (kFALSE),
  fPublishOutput (kTRUE),
  fProcessNJets (2),
  fJetsName (Names::gkPFJetBrn),
  fJetsFromBranch (kTRUE),
  fJets (0),
  fPfCandidatesName (Names::gkPFCandidatesBrn),
  fPfCandidatesFromBranch(kTRUE),
  fPfCandidates (0),
  fPileUpDenName(Names::gkPileupEnergyDensityBrn),
  fPileUpDenFromBranch(kTRUE),
  fPileUpDen(0),
  fVertexesName (ModNames::gkGoodVertexesName),
  fVertexesFromBranch(kFALSE),
  fVertexes(0),
  fXlFatJetsName ("XlFatJets"),
  fXlSubJetsName ("XlSubJets"),
  fSoftDropZCut (0.1),     
  fSoftDropMuCut (1.),  
  fPruneZCut (0.1),     
  fPruneDistCut (0.5),  
  fFilterN (3),      
  fFilterRad (0.2),     
  fTrimRad (0.05),       
  fTrimPtFrac (0.03),    
  fConeSize (0.6),
  fCounter (0)
{
  // Constructor.
}

FillerXlJets::~FillerXlJets()
{
  // Destructor
  if (fXlSubJets)
    delete fXlSubJets;
  if (fXlFatJets)
    delete fXlFatJets;
}

//--------------------------------------------------------------------------------------------------
void FillerXlJets::Process()
{
  // make sure the out collections are empty before starting
  fXlFatJets->Delete();  
  fXlSubJets->Delete();  
  
  // Load the branches we want to work with
  LoadEventObject(fJetsName,fJets,fJetsFromBranch);
  if (fQGTaggingActive) {
    LoadEventObject(fPileUpDenName,fPileUpDen,fPileUpDenFromBranch);
    LoadEventObject(fVertexesName,fVertexes,fVertexesFromBranch);
  }

  // Loop over PFCandidates and unmark them : necessary for skimming
  for (UInt_t i=0; i<fPfCandidates->GetEntries(); ++i) 
    fPfCandidates->At(i)->UnmarkMe();
 
  // Setup pileup density for QG computation
  if (fQGTaggingActive)
    fQGTagger->SetRhoIso(fPileUpDen->At(0)->RhoRandomLowEta());
  
  // Loop over jets
  for (UInt_t i=0; i<fJets->GetEntries(); ++i) {

    // consider only the first fProcessNJets jets
    if (i >= fProcessNJets)
      break; 
      
    const PFJet *jet = dynamic_cast<const PFJet*>(fJets->At(i));
    if (! jet) {
      printf(" FillerXlJets::Process() - ERROR - jets provided are not PFJets.");
      break;
    }
 
    // mark jet (and consequently its consituents) for further use in skim
    jet->Mark();       
    
    // perform Nsubjettiness analysis and fill the extended XlFatJet object
    // this method will also fill the SubJet collection       
    FillXlFatJet(jet);      
    
  }    
  
  return;
}

//--------------------------------------------------------------------------------------------------
void FillerXlJets::SlaveBegin()
{
  // Run startup code on the computer (slave) doing the actual analysis. 
  ReqEventObject(fJetsName,fJets,fJetsFromBranch);
  ReqEventObject(fPfCandidatesName,fPfCandidates,fPfCandidatesFromBranch);
  ReqEventObject(fPileUpDenName,fPileUpDen,fPileUpDenFromBranch);
  ReqEventObject(fVertexesName,fVertexes,fVertexesFromBranch);

  // Initialize area caculation (done with ghost particles)

  // Create the new output collection
  fXlFatJets = new XlFatJetArr(16,fXlFatJetsName);
  fXlSubJets = new XlSubJetArr(16,fXlSubJetsName);
  // Publish collection for further usage in the analysis
  if (fPublishOutput) {
    PublishObj(fXlFatJets);
    PublishObj(fXlSubJets);
  }

  // Prepare pruner
  fPruner = new fastjet::Pruner(fastjet::cambridge_algorithm,fPruneZCut,fPruneDistCut);
  // Prepare filterer
  fFilterer = new fastjet::Filter(fastjet::JetDefinition(fastjet::cambridge_algorithm,fFilterRad), 
                                  fastjet::SelectorNHardest(fFilterN));
  // Prepare trimmer
  fTrimmer = new fastjet::Filter(fastjet::Filter(fastjet::JetDefinition(fastjet::kt_algorithm,fTrimRad),
                                 fastjet::SelectorPtFractionMin(fTrimPtFrac)));
    
  // CA constructor (fConeSize = 0.6 for antiKt) - reproducing paper 1: http://arxiv.org/abs/1011.2268
  fCAJetDef = new fastjet::JetDefinition(fastjet::antikt_algorithm, fConeSize);
  
  // Initialize area caculation (done with ghost particles)
  int activeAreaRepeats = 1;
  double ghostArea = 0.01;
  double ghostEtaMax = 7.0;
  fActiveArea = new fastjet::GhostedAreaSpec(ghostEtaMax,activeAreaRepeats,ghostArea);
  fAreaDefinition = new fastjet::AreaDefinition(fastjet::active_area_explicit_ghosts,*fActiveArea);  
  
  // Initialize QGTagger class
  fQGTagger = new QGTagger(fQGTaggerCHS);
  
  return;
}

//--------------------------------------------------------------------------------------------------
void FillerXlJets::SlaveTerminate()
{
}

//--------------------------------------------------------------------------------------------------
void FillerXlJets::FillXlFatJet(const PFJet *pPFJet)
{
  // Prepare and store in an array a new FatJet 
  XlFatJet *fatJet = fXlFatJets->Allocate();
  new (fatJet) XlFatJet(*pPFJet);

  // Compute and store weighted charge
  fatJet->SetCharge();
  
  // Prepare and store QG tagging info
  float qgValue = -1.;
  if (fQGTaggingActive) {
    fQGTagger->CalculateVariables(pPFJet, fVertexes);
    qgValue = fQGTagger->QGValue();
  }
  fatJet->SetQGTag(qgValue);
    
  std::vector<fastjet::PseudoJet> fjParts;
  // Push all particle flow candidates of the input PFjet into fastjet particle collection
  for (UInt_t j=0; j<pPFJet->NPFCands(); ++j) {
    const PFCandidate *pfCand = pPFJet->PFCand(j);
    fjParts.push_back(fastjet::PseudoJet(pfCand->Px(),pfCand->Py(),pfCand->Pz(),pfCand->E()));
    fjParts.back().set_user_index(j);
  }	

  // Setup the cluster for fastjet
  fastjet::ClusterSequenceArea fjClustering(fjParts,*fCAJetDef,*fAreaDefinition);

  // ---- Fastjet is ready ----

  // Produce a new set of jets based on the fastjet particle collection and the defined clustering
  // Cut off fat jets with pt < 10 GeV and consider only the hardest jet of the output collection
  std::vector<fastjet::PseudoJet> fjOutJets = sorted_by_pt(fjClustering.inclusive_jets(10.)); 

  // Check that the output collection size is non-null, otherwise nothing to be done further
  if (fjOutJets.size() < 1) {
    printf(" FillerXlJets::FillXlFatJet() - WARNING - input PFJet produces null reclustering output");
    return;
  }
  fastjet::PseudoJet fjJet = fjOutJets[0];
    
  // Compute the subjettiness
  fastjet::contrib::Njettiness::AxesMode axisMode = fastjet::contrib::Njettiness::onepass_wta_kt_axes;
  fastjet::contrib::Njettiness::MeasureMode measureMode = fastjet::contrib::Njettiness::unnormalized_measure;
  double beta = 1.0;
  fastjet::contrib::Nsubjettiness  nSub1(1,axisMode,measureMode,beta);
  fastjet::contrib::Nsubjettiness  nSub2(2,axisMode,measureMode,beta);
  fastjet::contrib::Nsubjettiness  nSub3(3,axisMode,measureMode,beta);
  double tau1 = nSub1(fjJet);
  double tau2 = nSub2(fjJet);
  double tau3 = nSub3(fjJet);

  // Compute the energy correlation function ratios
  fastjet::contrib::EnergyCorrelatorDoubleRatio ECR2b0  (2,0. ,fastjet::contrib::EnergyCorrelator::pt_R);
  fastjet::contrib::EnergyCorrelatorDoubleRatio ECR2b0p2(2,0.2,fastjet::contrib::EnergyCorrelator::pt_R);
  fastjet::contrib::EnergyCorrelatorDoubleRatio ECR2b0p5(2,0.5,fastjet::contrib::EnergyCorrelator::pt_R);
  fastjet::contrib::EnergyCorrelatorDoubleRatio ECR2b1  (2,1.0,fastjet::contrib::EnergyCorrelator::pt_R);
  fastjet::contrib::EnergyCorrelatorDoubleRatio ECR2b2  (2,2.0,fastjet::contrib::EnergyCorrelator::pt_R);
  double C2b0   = ECR2b0(fjJet);
  double C2b0p2 = ECR2b0p2(fjJet);
  double C2b0p5 = ECR2b0p5(fjJet);
  double C2b1   = ECR2b1(fjJet);
  double C2b2   = ECR2b2(fjJet);

  // Compute Q-jets volatility
  std::vector<fastjet::PseudoJet> constits;
  GetJetConstituents(fjJet, constits, 0.01);
  double QJetVol = GetQjetVolatility(constits, 25, fCounter*25);
  fCounter++;
  constits.clear();

  // Compute groomed masses
  fastjet::contrib::SoftDropTagger softDropSDb0(0.0, fSoftDropZCut, fSoftDropMuCut);
  fastjet::contrib::SoftDropTagger softDropSDb2(2.0, fSoftDropZCut, fSoftDropMuCut);
  fastjet::contrib::SoftDropTagger softDropSDbm1(-1.0, fSoftDropZCut, fSoftDropMuCut);
  double MassSDb0 = (softDropSDb0(fjJet)).m();
  double MassSDb2 = (softDropSDb2(fjJet)).m();
  double MassSDbm1 = (softDropSDbm1(fjJet)).m();

  fastjet::PseudoJet fjJetPruned = (*fPruner)(fjJet);
  double MassPruned = fjJetPruned.m();
  double MassFiltered = ((*fFilterer)(fjJet)).m();
  double MassTrimmed = ((*fTrimmer)(fjJet)).m();
    
  // ---- Fastjet is done ----
          
  // Store the subjettiness values
  fatJet->SetTau1(tau1);
  fatJet->SetTau2(tau2);
  fatJet->SetTau3(tau3);

  // Store the energy correlation values
  fatJet->SetC2b0(C2b0);  
  fatJet->SetC2b0p2(C2b0p2);
  fatJet->SetC2b0p5(C2b0p5);
  fatJet->SetC2b1(C2b1);  
  fatJet->SetC2b2(C2b2);  

  // Store the Qjets volatility
  fatJet->SetQJetVol(QJetVol);  
  
  // Store the groomed masses
  fatJet->SetMassSDb0(MassSDb0);     
  fatJet->SetMassSDb2(MassSDb2);     
  fatJet->SetMassSDbm1(MassSDbm1);    
  fatJet->SetMassPruned(MassPruned);   
  fatJet->SetMassFiltered(MassFiltered);  
  fatJet->SetMassTrimmed(MassTrimmed);  

  // Store the color pull
  fatJet->SetPull(GetPull(fjJet,0.01).Mod());   
 
  // Loop on the subjets and fill the subjet Xl collections - do it according to the user request
  // In case of CA declustering subjet axis = subjet 
  if (fFillVSubJets) {
    std::vector<fastjet::PseudoJet> fjVSubJets;
    std::vector<fastjet::PseudoJet> fjVSubAxes;
    if (fNSubDeclustering) {
      fjVSubJets = nSub2.currentSubjets();
      fjVSubAxes = nSub2.currentAxes();
    }
    else {
      int nSubJPruned = std::min<unsigned int>(fjJetPruned.constituents().size(),2);
      fjVSubJets = fjJetPruned.associated_cluster_sequence()->exclusive_subjets(fjJetPruned,nSubJPruned);
      fjVSubAxes = fjVSubJets;    
    }
    std::vector<fastjet::PseudoJet> fjSubJetsSorted = sorted_by_pt(fjVSubJets);    
    // Store the color pull angle: either choose 2-prong or 3-prong subclustering!
    fatJet->SetPullAngle(GetPullAngle(fjSubJetsSorted,0.01));   
    FillXlSubJets(fjVSubJets,fjVSubAxes,fatJet,XlSubJet::ESubJetType::eV);
  } 
  if (fFillTopSubJets) {
    std::vector<fastjet::PseudoJet> fjTopSubJets;
    std::vector<fastjet::PseudoJet> fjTopSubAxes;
    if (fNSubDeclustering) {
      fjTopSubJets = nSub3.currentSubjets();
      fjTopSubAxes = nSub3.currentAxes();
    }
    else {
      int nSubJPruned = std::min<unsigned int>(fjJetPruned.constituents().size(),3);
      fjTopSubJets = fjJetPruned.associated_cluster_sequence()->exclusive_subjets(fjJetPruned,nSubJPruned);
      fjTopSubAxes = fjTopSubJets;    
    }
    std::vector<fastjet::PseudoJet> fjSubJetsSorted = sorted_by_pt(fjTopSubJets);    
    // Store the color pull angle: either choose 2-prong or 3-prong subclustering!
    fatJet->SetPullAngle(GetPullAngle(fjSubJetsSorted,0.01));   
    FillXlSubJets(fjTopSubJets,fjTopSubAxes,fatJet,XlSubJet::ESubJetType::eTop);
  } 
  // Sort subjets according to pt and add the subjets to the fat jet
  fXlSubJets->Sort();
  for (UInt_t iSubJ=0; iSubJ<fXlSubJets->GetEntries(); ++iSubJ)
    fatJet->AddSubJet(fXlSubJets->At(iSubJ));    
   
  return;
}

//--------------------------------------------------------------------------------------------------
void FillerXlJets::FillXlSubJets(std::vector<fastjet::PseudoJet> &fjSubJets, std::vector<fastjet::PseudoJet> &fjSubAxes,
                                 XlFatJet *pFatJet, XlSubJet::ESubJetType subJetType)
{
  for (int iSJet=0; iSJet < (int) fjSubJets.size(); iSJet++) {
    XlSubJet *subJet = fXlSubJets->Allocate();
    // Prepare and store in an array a new SubJet 
    new (subJet) XlSubJet(fjSubJets[iSJet].px(),
                          fjSubJets[iSJet].py(),
                          fjSubJets[iSJet].pz(),
                          fjSubJets[iSJet].e());

    // Store the subjet axis 
    ThreeVector subAxis(fjSubAxes[iSJet].px(),fjSubAxes[iSJet].py(),fjSubAxes[iSJet].pz());
    subJet->SetAxis(subAxis);

    // Store the QG tagging variable
    if (fQGTaggingActive)
      FillSubjetQGTagging(fjSubJets[iSJet], 0.01, subJet, pFatJet);
    
    // Store the subjet type value 
    subJet->SetSubJetType(subJetType);
                          
  }
    
  return;    
}

//--------------------------------------------------------------------------------------------------
void FillerXlJets::GetJetConstituents(fastjet::PseudoJet &jet, std::vector <fastjet::PseudoJet> &constits,  
                                      float constitsPtMin)
{
  // Loop on input jet constituents vector and discard very soft particles (ghosts)
  for (unsigned int iPart = 0; iPart < jet.constituents().size(); iPart++) {
    if (jet.constituents()[iPart].perp() < constitsPtMin)
      continue;
    constits.push_back(jet.constituents()[iPart]);        
  }
  
  return;
}

//--------------------------------------------------------------------------------------------------
double FillerXlJets::GetQjetVolatility(std::vector <fastjet::PseudoJet> &constits, int QJetsN, int seed)
{  
  std::vector<float> qjetmasses;
  
  double zcut(0.1), dcut_fctr(0.5), exp_min(0.), exp_max(0.), rigidity(0.1), truncationFactor(0.01);
  
  QjetsPlugin qjet_plugin(zcut, dcut_fctr, exp_min, exp_max, rigidity, truncationFactor);
  fastjet::JetDefinition qjet_def(&qjet_plugin);
  
  for(unsigned int ii = 0 ; ii < (unsigned int) QJetsN ; ii++){    
    qjet_plugin.SetRandSeed(seed+ii); // new feature in Qjets to set the random seed
    fastjet::ClusterSequence qjet_seq(constits, qjet_def);
    
    vector<fastjet::PseudoJet> inclusive_jets2 = sorted_by_pt(qjet_seq.inclusive_jets(5.0));
    if (inclusive_jets2.size()>0) { qjetmasses.push_back( inclusive_jets2[0].m() ); }          
  }
  
  // find RMS of a vector
  float qjetsRMS = FindRMS( qjetmasses );
  // find mean of a vector
  float qjetsMean = FindMean( qjetmasses );
  float qjetsVolatility = qjetsRMS/qjetsMean;
  return qjetsVolatility;
}

//--------------------------------------------------------------------------------------------------
double FillerXlJets::FindRMS(std::vector<float> qjetmasses)
{
  float total = 0.;
  float ctr = 0.;
  for (unsigned int i = 0; i < qjetmasses.size(); i++){
      total = total + qjetmasses[i];
      ctr++;
  }
  float mean = total/ctr;
  
  float totalsquared = 0.;
  for (unsigned int i = 0; i < qjetmasses.size(); i++){
    totalsquared += (qjetmasses[i] - mean)*(qjetmasses[i] - mean) ;
  }
  float RMS = sqrt( totalsquared/ctr );
  return RMS;
}

//--------------------------------------------------------------------------------------------------
double FillerXlJets::FindMean(std::vector<float> qjetmasses)
{
  float total = 0.;
  float ctr = 0.;
  for (unsigned int i = 0; i < qjetmasses.size(); i++){
      total = total + qjetmasses[i];
      ctr++;
  }
  return total/ctr;
}

//--------------------------------------------------------------------------------------------------
void FillerXlJets::FillSubjetQGTagging(fastjet::PseudoJet &jet, float constitsPtMin, 
                                       XlSubJet *pSubJet, XlFatJet *pFatJet)
{
  // Prepare a PFJet to compute the QGTagging
  PFJet pfJet(jet.px(),jet.py(),jet.pz(),jet.e());

  // Loop on input jet constituents vector and discard very soft particles (ghosts)
  for (unsigned int iPart = 0; iPart < jet.constituents().size(); iPart++) {
    if (jet.constituents()[iPart].perp() < constitsPtMin)
      continue;
    int thisPFCandIndex = jet.constituents()[iPart].user_index();
    // Add the constituent to the PF subjet
    pfJet.AddPFCand(pFatJet->PFCand(thisPFCandIndex));      
  }

  // Compute the subjet QGTagging
  if (fQGTaggingActive) {
    fQGTagger->CalculateVariables(&pfJet, fVertexes);
    pSubJet->SetQGTag(fQGTagger->QGValue());
    pSubJet->SetQGPtD(fQGTagger->GetPtD());
    pSubJet->SetQGAxis1(fQGTagger->GetAxis1());
    pSubJet->SetQGAxis2(fQGTagger->GetAxis2());
    pSubJet->SetQGMult(fQGTagger->GetMult());
  }

  return;
}

//--------------------------------------------------------------------------------------------------
TVector2 FillerXlJets::GetPull(fastjet::PseudoJet &jet, float constitsPtMin)
{
  double dYSum   = 0;
  double dPhiSum = 0;
  // Loop on input jet constituents vector and discard very soft particles (ghosts)
  for (unsigned int iPart = 0; iPart < jet.constituents().size(); iPart++) {
    if (jet.constituents()[iPart].perp() < constitsPtMin)
      continue;
    double dY     = jet.constituents()[iPart].rapidity()-jet.rapidity();
    double dPhi   = MathUtils::DeltaPhi(jet.constituents()[iPart].phi(),jet.phi());
    double weight = jet.constituents()[iPart].pt()*sqrt(dY*dY + dPhi*dPhi);
    dYSum   += weight*dY;
    dPhiSum += weight*dPhi;
  }
  return TVector2(dYSum/jet.pt(), dPhiSum/jet.pt());
}

//--------------------------------------------------------------------------------------------------
double FillerXlJets::GetPullAngle(std::vector<fastjet::PseudoJet> &fjSubJets, float constitsPtMin)
{
  // Subject collection already sorted by pt
  // Consider only the leading and the subleading for the pull angle computation
  // work in dy-dphi space of leading subjet

  // Exclude cases where there is no second subjet (input jet made by one particle)
  if (fjSubJets.size() < 2)
    return -20.;

  TVector2 lPull = GetPull(fjSubJets[0],constitsPtMin);
  TVector2 lJet(fjSubJets[1].rapidity()-fjSubJets[0].rapidity(), 
                MathUtils::DeltaPhi(fjSubJets[1].phi(), fjSubJets[0].phi()));
  double lThetaP = lPull.DeltaPhi(lJet);
  return lThetaP;
}

#include "MitPhysics/Init/interface/ModNames.h"
#include "MitMonoJet/TreeFiller/interface/FillerXlJets.h"
#include "MitAna/DataTree/interface/PFJetCol.h"

#include "MitMonoJet/DataTree/interface/XlJet.h"
#include "MitCommon/DataFormats/interface/Types.h"
#include "MitCommon/MathTools/interface/MathUtils.h"

using namespace mithep;

ClassImp(mithep::FillerXlJets)

//--------------------------------------------------------------------------------------------------
FillerXlJets::FillerXlJets(const char *name, const char *title) :
  BaseMod (name,title),
  fIsData (kTRUE),
  fQGTaggingActive (kTRUE),
  fQGTaggerCHS (kFALSE),
  fPublishOutput (kTRUE),
  fJetsName (Names::gkPFJetBrn),
  fJetsFromBranch (kTRUE),
  fJets (0),
  fPileUpDenName(Names::gkPileupEnergyDensityBrn),
  fPileUpDenFromBranch(kTRUE),
  fPileUpDen(0),
  fVertexesName (ModNames::gkGoodVertexesName),
  fVertexesFromBranch(kFALSE),
  fVertexes(0),
  fXlJetsName ("XlJets")
{
  // Constructor.
}

FillerXlJets::~FillerXlJets()
{
  // Destructor
  if (fXlJets)
    delete fXlJets;
  
  delete fQGTagger;
}

//--------------------------------------------------------------------------------------------------
void FillerXlJets::Process()
{
  // make sure the out collections are empty before starting
  fXlJets->Delete();  
  
  // Load the branches we want to work with
  LoadEventObject(fJetsName,fJets,fJetsFromBranch);
  if (fQGTaggingActive) {
    LoadEventObject(fPileUpDenName,fPileUpDen,fPileUpDenFromBranch);
    LoadEventObject(fVertexesName,fVertexes,fVertexesFromBranch);
  }
 
  // Setup pileup density for QG computation
  if (fQGTaggingActive)
    fQGTagger->SetRhoIso(fPileUpDen->At(0)->RhoRandomLowEta());
  
  // Loop over jets
  for (UInt_t i=0; i<fJets->GetEntries(); ++i) {
      
    const PFJet *jet = dynamic_cast<const PFJet*>(fJets->At(i));
    if (! jet) {
      printf(" FillerXlJets::Process() - ERROR - jets provided are not PFJets.");
      break;
    }
 
    // mark jet (and consequently its consituents) for further use in skim
    jet->Mark();       
    
    // perform jet analysis and fill the extended XlJet object
    FillXlJet(jet);      
    
  }    
  // Trim the output collection
  fXlJets->Trim();
  
  return;
}

//--------------------------------------------------------------------------------------------------
void FillerXlJets::SlaveBegin()
{
  // Run startup code on the computer (slave) doing the actual analysis. 
  ReqEventObject(fJetsName,fJets,fJetsFromBranch);
  ReqEventObject(fPileUpDenName,fPileUpDen,fPileUpDenFromBranch);
  ReqEventObject(fVertexesName,fVertexes,fVertexesFromBranch);

  // Initialize area caculation (done with ghost particles)

  // Create the new output collection
  fXlJets = new XlJetArr(16,fXlJetsName);
  // Publish collection for further usage in the analysis
  if (fPublishOutput)
    PublishObj(fXlJets);
  
  // Initialize QGTagger class
  fQGTagger = new QGTagger(fQGTaggerCHS);
  
  return;
}

//--------------------------------------------------------------------------------------------------
void FillerXlJets::SlaveTerminate()
{
}

//--------------------------------------------------------------------------------------------------
void FillerXlJets::FillXlJet(const PFJet *pPFJet)
{
  // Prepare and store in an array a new FatJet 
  XlJet *newXlJet = fXlJets->Allocate();
  new (newXlJet) XlJet(*pPFJet);

  // Compute and store weighted charge
  newXlJet->SetCharge();
  
  // Prepare and store QG tagging info
  float qgValue = -1.;
  if (fQGTaggingActive) {
    fQGTagger->CalculateVariables(pPFJet, fVertexes);
    qgValue = fQGTagger->QGValue();
  }
<<<<<<< HEAD
  fatJet->SetQGTag(qgValue);
    
  std::vector<fastjet::PseudoJet> fjParts;
  // Push all particle flow candidates of the input PFjet into fastjet particle collection
  for (UInt_t j=0; j<pPFJet->NPFCands(); ++j) {
    const PFCandidate *pfCand = pPFJet->PFCand(j);
    fjParts.push_back(fastjet::PseudoJet(pfCand->Px(),pfCand->Py(),pfCand->Pz(),pfCand->E()));
    fjParts.back().set_user_index(j);
  }	

  // Setup the cluster for fastjet
  fastjet::ClusterSequenceArea *fjClustering =
    new fastjet::ClusterSequenceArea(fjParts,*fCAJetDef,*fAreaDefinition);

  // ---- Fastjet is ready ----

  // Produce a new set of jets based on the fastjet particle collection and the defined clustering
  // Cut off fat jets with pt < 10 GeV and consider only the hardest jet of the output collection
  std::vector<fastjet::PseudoJet> fjOutJets;
  if (fFillTopSubJets) fjOutJets = sorted_by_pt(fjClustering->inclusive_jets(0.)); 
  else fjOutJets = sorted_by_pt(fjClustering->inclusive_jets(10.)); 

  // Check that the output collection size is non-null, otherwise nothing to be done further
  if (fjOutJets.size() < 1) {
    printf(" FillerXlJets::FillXlFatJet() - WARNING - input PFJet produces null reclustering output. skipping event!");

    fjClustering->delete_self_when_unused();
    delete fjClustering;

    this->SkipEvent(); 
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
  fastjet::contrib::SoftDropTagger softDropSDb1(1.0, fSoftDropZCut, fSoftDropMuCut);
  fastjet::contrib::SoftDropTagger softDropSDb2(2.0, fSoftDropZCut, fSoftDropMuCut);
  fastjet::contrib::SoftDropTagger softDropSDbm1(-1.0, fSoftDropZCut, fSoftDropMuCut);
  double MassSDb0 = (softDropSDb0(fjJet)).m();
  double MassSDb1 = (softDropSDb1(fjJet)).m();
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
  fatJet->SetMassSDb1(MassSDb1);     
  fatJet->SetMassSDb2(MassSDb2);     
  fatJet->SetMassSDbm1(MassSDbm1);    
  fatJet->SetMassPruned(MassPruned);   
  fatJet->SetMassFiltered(MassFiltered);  
  fatJet->SetMassTrimmed(MassTrimmed);  

  // Store the color pull
  fatJet->SetPull(GetPull(fjJet,0.01).Mod());   
 
  // Loop on the subjets and fill the subjet Xl collections - do it according to the user request
  if (fFillVSubJets) {
    std::vector<fastjet::PseudoJet> fjVSubJets;
    if (fNSubDeclustering)
      fjVSubJets = nSub2.currentSubjets();
    else {
      int nSubJPruned = std::min<unsigned int>(fjJetPruned.constituents().size(),2);
      fjVSubJets = fjJetPruned.associated_cluster_sequence()->exclusive_subjets(fjJetPruned,nSubJPruned);
    }
    // Order the subjets according to their pt and discard zero pt subjets
    std::vector<fastjet::PseudoJet> fjSubJetsSorted = Sorted_by_pt_min_pt(fjVSubJets,0.01);    
    // Store the color pull angle: either choose 2-prong or 3-prong subclustering!
    fatJet->SetPullAngle(GetPullAngle(fjSubJetsSorted,0.01));   
    FillXlSubJets(fjSubJetsSorted,fatJet,XlSubJet::ESubJetType::eV);
  } 
  if (fFillTopSubJets) {
    std::vector<fastjet::PseudoJet> fjTopSubJets;
    if (fNSubDeclustering) 
      fjTopSubJets = nSub3.currentSubjets();
    else {
      int nSubJPruned = std::min<unsigned int>(fjJetPruned.constituents().size(),3);
      fjTopSubJets = fjJetPruned.associated_cluster_sequence()->exclusive_subjets(fjJetPruned,nSubJPruned);
    }
    // Order the subjets according to their pt
    std::vector<fastjet::PseudoJet> fjSubJetsSorted = Sorted_by_pt_min_pt(fjTopSubJets,0.01);    
    // Store the color pull angle: either choose 2-prong or 3-prong subclustering!
    fatJet->SetPullAngle(GetPullAngle(fjSubJetsSorted,0.01));   
    FillXlSubJets(fjSubJetsSorted,fatJet,XlSubJet::ESubJetType::eTop);
  } 
  // Trim the output collections
  fXlSubJets->Trim();
  fXlFatJets->Trim();
   
  // Memory cleanup
  fjClustering->delete_self_when_unused();
  delete fjClustering;
   
  return;
}

//--------------------------------------------------------------------------------------------------
void FillerXlJets::FillXlSubJets(std::vector<fastjet::PseudoJet> &fjSubJets,
                                 XlFatJet *pFatJet, XlSubJet::ESubJetType subJetType)
{
  for (int iSJet=0; iSJet < (int) fjSubJets.size(); iSJet++) {
    XlSubJet *subJet = fXlSubJets->Allocate();
    // Prepare and store in an array a new SubJet 
    new (subJet) XlSubJet(fjSubJets[iSJet].px(),
                          fjSubJets[iSJet].py(),
                          fjSubJets[iSJet].pz(),
                          fjSubJets[iSJet].e());

    // Store the QG tagging variable
    if (fQGTaggingActive)
      FillSubjetQGTagging(fjSubJets[iSJet], 0.01, subJet, pFatJet);
    
    // Store the subjet type value 
    subJet->SetSubJetType(subJetType);

    // Add the subjet to the relative fatjet 
    pFatJet->AddSubJet(subJet);                              

  }
=======
  newXlJet->SetQGTag(qgValue);
>>>>>>> 36adccb0a04e529411a9763ba098d7bc8feffc3f
    
  // Prepare and store jet pull info
  TVector2 newXlJetPull = GetPull(pPFJet);
  newXlJet->SetPullY(newXlJetPull.X());
  newXlJet->SetPullPhi(newXlJetPull.Y());
      
  return;
}

//--------------------------------------------------------------------------------------------------
TVector2 FillerXlJets::GetPull(const PFJet *inPFJet)
{
  double dYSum   = 0;
  double dPhiSum = 0;
  const unsigned int nPFCands = inPFJet->NPFCands();

  // Loop on input jet constituents and get the color pull  
  for(unsigned int ipf=0; ipf<nPFCands; ipf++) {    
    const PFCandidate *pfCand = inPFJet->PFCand(ipf);
    double pt_i=0, y_i=0, phi_i=0;
    pt_i = pfCand->Pt();
    y_i = pfCand->Rapidity();
    phi_i = pfCand->Phi();

    double dY   = y_i - inPFJet->Rapidity();
    double dPhi = MathUtils::DeltaPhi(phi_i,inPFJet->Phi());
    double weight = pt_i*sqrt(dY*dY + dPhi*dPhi);
    dYSum += weight*dY;
    dPhiSum += weight*dPhi;
  }

  return TVector2(dYSum/inPFJet->Pt(), dPhiSum/inPFJet->Pt());    
}

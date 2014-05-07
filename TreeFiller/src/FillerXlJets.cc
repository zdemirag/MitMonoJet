#include "MitPhysics/Init/interface/ModNames.h"
#include "MitMonoJet/TreeFiller/interface/FillerXlJets.h"
#include "MitAna/DataTree/interface/PFJetCol.h"

#include "MitMonoJet/DataTree/interface/XlSubJet.h"
#include "MitMonoJet/DataTree/interface/XlFatJet.h"
#include "MitCommon/DataFormats/interface/Vect4M.h"
#include "MitCommon/DataFormats/interface/Vect3.h"
#include "MitCommon/DataFormats/interface/Types.h"

using namespace mithep;

ClassImp(mithep::FillerXlJets)

//--------------------------------------------------------------------------------------------------
FillerXlJets::FillerXlJets(const char *name, const char *title) :
  BaseMod (name,title),
  fIsData (kTRUE),
  fFillVSubJets (kTRUE),
  fFillTopSubJets (kFALSE),
  fBTaggingActive (kFALSE),
  fQGTaggingActive (kFALSE),
  fPublishOutput (kTRUE),
  fProcessNJets (2),
  fJetsName (Names::gkPFJetBrn),
  fJetsFromBranch (kTRUE),
  fJets (0),
  fPfCandidatesName (Names::gkPFCandidatesBrn),
  fPfCandidatesFromBranch(kTRUE),
  fPfCandidates (0),
  fXlFatJetsName ("XlFatJets"),
  fXlSubJetsName ("XlSubJets"),
  fPrune (kFALSE),         
  fFilter (kFALSE),        
  fTrim (kFALSE),          
  fPruneZCut (0.1),     
  fPruneDistCut (0.5),  
  fFilterN (3),      
  fFilterRad (0.2),     
  fTrimRad (0.05),       
  fTrimPtFrac (0.03),    
  fConeSize (0.6)
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

  // Loop over PFCandidates and unmark them : necessary for skimming
  for (UInt_t i=0; i<fPfCandidates->GetEntries(); ++i) 
    fPfCandidates->At(i)->UnmarkMe();
 
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
  // Run startup code on the computer (slave) doing the actual analysis. Here, we just request the
  // particle flow collection branch.
  ReqEventObject(fJetsName,fJets,fJetsFromBranch);
  ReqEventObject(fPfCandidatesName,fPfCandidates,fPfCandidatesFromBranch);

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
  
  return;
}

//--------------------------------------------------------------------------------------------------
void FillerXlJets::SlaveTerminate()
{
}

//--------------------------------------------------------------------------------------------------
void FillerXlJets::FillXlFatJet(const PFJet *pPFJet)
{
  
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

  // If required by the user apply grooming algorithm on the jet
  fastjet::PseudoJet fjGroomedJet;
  if (fPrune)
    fjGroomedJet = (*fPruner)(fjOutJets[0]);
  else if (fFilter)
    fjGroomedJet = (*fFilterer)(fjOutJets[0]);
  else if (fTrim)
    fjGroomedJet = (*fTrimmer)(fjOutJets[0]);
  else 
    fjGroomedJet = fjOutJets[0];

  // Compute the subjettiness
  fastjet::contrib::Njettiness::AxesMode axisMode = fastjet::contrib::Njettiness::onepass_wta_kt_axes;
  fastjet::contrib::Njettiness::MeasureMode measureMode = fastjet::contrib::Njettiness::unnormalized_measure;
  double beta = 1.0;
  fastjet::contrib::Nsubjettiness  nSub1(1,axisMode,measureMode,beta);
  fastjet::contrib::Nsubjettiness  nSub2(2,axisMode,measureMode,beta);
  fastjet::contrib::Nsubjettiness  nSub3(3,axisMode,measureMode,beta);
  double tau1 = nSub1(fjGroomedJet);
  double tau2 = nSub2(fjGroomedJet);
  double tau3 = nSub3(fjGroomedJet);

  // ---- Fastjet is done ----
      
  // Prepare and store in an array a new FatJet 
  XlFatJet *fatJet = fXlFatJets->Allocate();
  new (fatJet) XlFatJet(*pPFJet);

  // Store the groomed momentum
  FourVectorM groomedMom(fjGroomedJet.pt(),fjGroomedJet.eta(),fjGroomedJet.phi(),fjGroomedJet.m());
  fatJet->SetGroomedMom(groomedMom);
    
  // Store the subjettiness values
  fatJet->SetTau1(tau1);
  fatJet->SetTau2(tau2);
  fatJet->SetTau3(tau3);

  // Loop on the subjets and fill the subjet Xl collections - do it according to the user request
  if (fFillVSubJets) {
    std::vector<fastjet::PseudoJet> fjVSubJets = nSub2.currentSubjets();
    std::vector<fastjet::PseudoJet> fjVSubAxes = nSub2.currentAxes();
    FillXlSubJets(fjVSubJets,fjVSubAxes,fatJet,XlSubJet::ESubJetType::eV);
  } // End scope of V-subjets filling
  if (fFillTopSubJets) {
    std::vector<fastjet::PseudoJet> fjTopSubJets = nSub3.currentSubjets();
    std::vector<fastjet::PseudoJet> fjTopSubAxes = nSub3.currentAxes();
    FillXlSubJets(fjTopSubJets,fjTopSubAxes,fatJet,XlSubJet::ESubJetType::eTop);
  } // End scope of Top-subjets filling
  
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
    
    // Store the subjet type value 
    subJet->SetSubJetType(subJetType);
                          
    // Add the subjet to the fatjet
    pFatJet->AddSubJet(subJet);
  }
    
  return;    
}

#include "MitPhysics/Init/interface/ModNames.h"
#include "MitMonoJet/TreeFiller/interface/FillerXlJets.h"
#include "MitAna/DataTree/interface/PFJetCol.h"

#include "MitMonoJet/DataTree/interface/XlSubJet.h"
#include "MitMonoJet/DataTree/interface/XlFatJet.h"
#include "MitCommon/DataFormats/interface/Vect4M.h"
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
 
  // Loop over jets and perform Nsubjettiness analysis (for now just stick with the first two jets)
  std::vector<fastjet::PseudoJet> lFjParts;
  for (UInt_t i=0; i<fJets->GetEntries(); ++i) {

    // consider only the first fProcessNJets jets
    if (i >= fProcessNJets)
      break; 
      
    const PFJet *jet = dynamic_cast<const PFJet*>(fJets->At(i));
    if (! jet) {
      printf(" FillerXlJets::Process() - ERROR - jets provided are not PFJets.");
      break;
    }
        
    // Push all particle flow candidates into fastjet particle collection
    for (UInt_t j=0; j<jet->NPFCands(); ++j) {
      const PFCandidate *pfCand = jet->PFCand(j);
      lFjParts.push_back(fastjet::PseudoJet(pfCand->Px(),pfCand->Py(),pfCand->Pz(),pfCand->E()));
      lFjParts.back().set_user_index(j);      
    }	
  }

  // ---- Fastjet is ready ----

  // Setup the cluster for fastjet
  fastjet::ClusterSequenceArea *lClustering =
    new fastjet::ClusterSequenceArea(lFjParts,*fCAJetDef,*fAreaDefinition);

  // Produce a new set of jets based on the fastjet particle collection and the defined clustering
  // Cut off fat jets with pt < 10 GeV
  std::vector<fastjet::PseudoJet> lOutJets = sorted_by_pt(lClustering->inclusive_jets(10.)); 
  
  // Fill the new collections with the output of fastjet
  FillfXlFatJets(lOutJets); // this method will also fill the SubJet collection
   
  // ---- Fastjet is done ----
      
  // Always cleanup
  if (lClustering)
    delete lClustering;
    
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
  
}

//--------------------------------------------------------------------------------------------------
void FillerXlJets::SlaveTerminate()
{
}

//--------------------------------------------------------------------------------------------------
void FillerXlJets::FillfXlFatJets(std::vector<fastjet::PseudoJet> &fjFatJets)
{
  for (int iJet=0; iJet < (int) fjFatJets.size(); iJet++) {
    // Skip very soft jets produced by fastjet clustering 
    if (fjFatJets[iJet].perp() < 10) continue;

    // If required by the user apply grooming algorithm on the jet
    fastjet::PseudoJet pJet;
    if (fPrune)
      pJet = (*fPruner)(fjFatJets[iJet]);
    else if (fFilter)
      pJet = (*fFilterer)(fjFatJets[iJet]);
    else if (fTrim)
      pJet = (*fTrimmer)(fjFatJets[iJet]);
    else 
      pJet = fjFatJets[iJet];
      
    // Prepare and store in an array a new FatJet 
    XlFatJet *fatJet = fXlFatJets->Allocate();
    new (fatJet) XlFatJet(pJet.px(),
                          pJet.py(),
                          pJet.pz(),
                          pJet.e());
    
    // Compute the subjettiness
    fastjet::contrib::Njettiness::AxesMode axisMode = fastjet::contrib::Njettiness::onepass_wta_kt_axes;
    fastjet::contrib::Njettiness::MeasureMode measureMode = fastjet::contrib::Njettiness::unnormalized_measure;
    double beta = 1.0;
    fastjet::contrib::Nsubjettiness  nSub1(1,axisMode,measureMode,beta);
    fastjet::contrib::Nsubjettiness  nSub2(2,axisMode,measureMode,beta);
    fastjet::contrib::Nsubjettiness  nSub3(3,axisMode,measureMode,beta);
    double tau1 = nSub1(pJet);
    double tau2 = nSub2(pJet);
    double tau3 = nSub3(pJet);

    // Store the subjettiness values
    fatJet->SetTau1(tau1);
    fatJet->SetTau2(tau2);
    fatJet->SetTau3(tau3);

    // Loop on the subjets and fill the subjet Xl collections - do it according to the user request
    if (fFillVSubJets) {
      std::vector<fastjet::PseudoJet> fjVSubJets = nSub2.currentSubjets();
      FillfXlSubJets(fjVSubJets,fatJet,XlSubJet::ESubJetType::eV);
    } // End scope of V-subjets filling
    if (fFillTopSubJets) {
      std::vector<fastjet::PseudoJet> fjTopSubJets = nSub3.currentSubjets();
      FillfXlSubJets(fjTopSubJets,fatJet,XlSubJet::ESubJetType::eTop);
    } // End scope of Top-subjets filling
  }
  
  return;
}

//--------------------------------------------------------------------------------------------------
void FillerXlJets::FillfXlSubJets(std::vector<fastjet::PseudoJet> &fjSubJets, XlFatJet *pFatJet,
                                 XlSubJet::ESubJetType subJetType)
{
  for (int iSJet=0; iSJet < (int) fjSubJets.size(); iSJet++) {
    XlSubJet *subJet = fXlSubJets->Allocate();
    // Prepare and store in an array a new SubJet 
    new (subJet) XlSubJet(fjSubJets[iSJet].px(),
                          fjSubJets[iSJet].py(),
                          fjSubJets[iSJet].pz(),
                          fjSubJets[iSJet].e());

    // Store the subjet type value 
    subJet->SetSubJetType(subJetType);
                          
    // Add the subjet to the fatjet
    pFatJet->AddSubJet(subJet);
  }
    
  return;    
}

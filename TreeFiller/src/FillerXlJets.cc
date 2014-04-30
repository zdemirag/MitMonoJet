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
  fBTaggingActive (kFALSE),
  fQGTaggingActive (kFALSE),
  fPublishOutput (kTRUE),
  fJetsName (Names::gkPFJetBrn),
  fJetsFromBranch (kTRUE),
  fJets (0),
  fPfCandidatesName (Names::gkPFCandidatesBrn),
  fPfCandidatesFromBranch(kTRUE),
  fPfCandidates (0),
  XlFatJetsName ("XlFatJets"),
  XlFatJets (new XlFatJetArr(16,XlFatJetsName)),
  XlSubJetsName ("XlSubJets"),
  XlSubJets (new XlSubJetArr(16,XlSubJetsName)),
  fConeSize (0.6),
  fPrune (1.)
{
  // Constructor.
}

FillerXlJets::~FillerXlJets()
{
  // Destructor
  if (XlSubJets)
    delete XlSubJets;
  if (XlFatJets)
    delete XlFatJets;
}

//--------------------------------------------------------------------------------------------------
void FillerXlJets::Process()
{
  // make sure the out collections are empty before starting
  XlFatJets->Delete();  
  XlSubJets->Delete();  
  
  // Load the branches we want to work with
  LoadEventObject(fJetsName,fJets,fJetsFromBranch);
 
  // Loop over jets and perform Nsubjettiness analysis (for now just stick with the first two jets)
  std::vector<fastjet::PseudoJet> lFjParts;
  for (UInt_t i=0; i<fJets->GetEntries(); ++i) {
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
    if (i > 0)
      break; // this is for now just considering the first two jets (consider all in the future)
  }

  // ---- Fastjet is ready ----

  // Setup the cluster for fastjet
  fastjet::ClusterSequenceArea *lClustering =
    new fastjet::ClusterSequenceArea(lFjParts,*fCAJetDef,*fAreaDefinition);

  // Produce a new set of jets based on the fastjet particle collection and the defined clustering
  std::vector<fastjet::PseudoJet> lOutJets = sorted_by_pt(lClustering->inclusive_jets());
  
  // Fill the new collections with the output of fastjet
  FillXlFatJets(lOutJets); // this method will also fill the SubJet collection
   
  // Access the 2 hardest jets as applicable
  fastjet::PseudoJet lJet1, lJet2;
  if (fPrune) {
    if (lOutJets.size() > 0)
      lJet1 = (*fPruner)(lOutJets[0]);
    if (lOutJets.size() > 1)
      lJet2 = (*fPruner)(lOutJets[1]);
  }
  else {
    if (lOutJets.size() > 0)
      lJet1 = lOutJets[0];
    if (lOutJets.size() > 1)
      lJet2 = lOutJets[1];
  }

  // ---- Fastjet is already done ----

  // Skip event if no jets are found
  if (lOutJets.size() == 0) {
    if (lClustering)
      delete lClustering;
    return;
  }
      
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

  // Default pruning parameters
  fPruner = new fastjet::Pruner(fastjet::cambridge_algorithm,0.1,0.5); // CMS Default
  
  // CA constructor (fConeSize = 0.6 for antiKt) - reproducing paper 1: http://arxiv.org/abs/1011.2268
  fCAJetDef = new fastjet::JetDefinition(fastjet::antikt_algorithm, fConeSize);
  
  // Initialize area caculation (done with ghost particles)
  int activeAreaRepeats = 1;
  double ghostArea = 0.01;
  double ghostEtaMax = 7.0;
  fActiveArea = new fastjet::GhostedAreaSpec(ghostEtaMax,activeAreaRepeats,ghostArea);
  fAreaDefinition = new fastjet::AreaDefinition(fastjet::active_area_explicit_ghosts,*fActiveArea);
  
  // Publish collection for further usage in the analysis
  if (fPublishOutput) {
    PublishObj(XlFatJets);
    PublishObj(XlSubJets);
  }

}

//--------------------------------------------------------------------------------------------------
void FillerXlJets::SlaveTerminate()
{
}

//--------------------------------------------------------------------------------------------------
void FillerXlJets::FillXlFatJets(std::vector<fastjet::PseudoJet> &fjFatJets)
{
  for (int iJet=0; iJet < (int) fjFatJets.size(); iJet++) {
    // Skip very soft jets produced by fastjet clustering 
    if (fjFatJets[iJet].perp() < 10) continue;

    // Prepare and store in an array a new FatJet 
    XlFatJet *fatJet = XlFatJets->Allocate();
    new (fatJet) XlFatJet(fjFatJets[iJet].px(),
                          fjFatJets[iJet].py(),
                          fjFatJets[iJet].pz(),
                          fjFatJets[iJet].e());
    
    // Compute the subjettiness
    fastjet::contrib::Njettiness::AxesMode axisMode = fastjet::contrib::Njettiness::onepass_wta_kt_axes;
    fastjet::contrib::Njettiness::MeasureMode measureMode = fastjet::contrib::Njettiness::unnormalized_measure;
    double beta = 1.0;
    fastjet::contrib::Nsubjettiness  nSub1(1,axisMode,measureMode,beta);
    fastjet::contrib::Nsubjettiness  nSub2(2,axisMode,measureMode,beta);
    fastjet::contrib::Nsubjettiness  nSub3(3,axisMode,measureMode,beta);
    double tau1 = nSub1(fjFatJets[iJet]);
    double tau2 = nSub2(fjFatJets[iJet]);
    double tau3 = nSub3(fjFatJets[iJet]);

    // Store the subjettiness values
    fatJet->SetTau1(tau1);
    fatJet->SetTau2(tau2);
    fatJet->SetTau3(tau3);

    // Loop on the subjets and fill the subjet Xl collection - for now only considering 2-prom structure
    std::vector<fastjet::PseudoJet> fjSubJets = nSub2.currentSubjets();
    FillXlSubJets(fjSubJets,fatJet);
  }
  
  return;
}

//--------------------------------------------------------------------------------------------------
void FillerXlJets::FillXlSubJets(std::vector<fastjet::PseudoJet> &fjSubJets, XlFatJet *pFatJet)
{
  for (int iSJet=0; iSJet < (int) fjSubJets.size(); iSJet++) {
    XlSubJet *subJet = XlSubJets->Allocate();
    // Prepare and store in an array a new SubJet 
    new (subJet) XlSubJet(fjSubJets[iSJet].px(),
                          fjSubJets[iSJet].py(),
                          fjSubJets[iSJet].pz(),
                          fjSubJets[iSJet].e());
                          
    // Add the subjet to the fatjet
    pFatJet->AddSubJet(subJet);
  }
    
  return;    
}

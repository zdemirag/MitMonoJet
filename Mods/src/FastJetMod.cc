#include "MitPhysics/Init/interface/ModNames.h"
#include "MitMonoJet/Mods/interface/FastJetMod.h"
#include "MitAna/DataTree/interface/PFJetCol.h"

#include "MitMonoJet/DataTree/interface/XlSubJet.h"
#include "MitMonoJet/DataTree/interface/XlFatJet.h"
#include "MitCommon/DataFormats/interface/Vect4M.h"
#include "MitCommon/DataFormats/interface/Vect3.h"
#include "MitCommon/DataFormats/interface/Types.h"
#include "MitCommon/MathTools/interface/MathUtils.h"

using namespace mithep;

ClassImp(mithep::FastJetMod)

//--------------------------------------------------------------------------------------------------
FastJetMod::FastJetMod(const char *name, const char *title) :
  BaseMod (name,title),
  fGetMatchBtag (kTRUE),
  fUseBambuJets (kFALSE),
  fBtaggedJetsName (Names::gkPFJetBrn),
  fBtaggedJetsFromBranch (kTRUE),
  fBtaggedJets (0),
  fJetsName (Names::gkPFJetBrn),
  fJetsFromBranch (kTRUE),
  fJets (0),
  fPfCandidatesName (Names::gkPFCandidatesBrn),
  fPfCandidatesFromBranch(kTRUE),
  fPfCandidates (0),
  fOutputJetsName ("AK8FJCHS"),
  fOutputJets (0),
  fJetConeSize (0.5),
  fFatJetConeSize (0.8),
  fParticleMinPt (0.001),
  fJetMinPt (20)
{
  // Constructor.
}

FastJetMod::~FastJetMod()
{
  // Destructor  
  delete fOutputJets;
  delete fJetDef;
  
  delete fActiveArea;
  delete fAreaDefinition;  
}

//--------------------------------------------------------------------------------------------------
void FastJetMod::Process()
{
  
  // Get input collections
  if (fGetMatchBtag) {
    LoadEventObject(fBtaggedJetsName,fBtaggedJets,fBtaggedJetsFromBranch);
    if (!fBtaggedJets) {
      SendError(kAbortModule,"Process","Pointer to input b-tagged jet collection %s null.",fBtaggedJetsName.Data());
      return;
    }
  }
  if (fUseBambuJets) {
    LoadEventObject(fJetsName,fJets,fJetsFromBranch);
    if (!fJets) {
      SendError(kAbortModule,"Process","Pointer to input jet collection %s null.",fJetsName.Data());
      return;
    }
  }
  LoadEventObject(fPfCandidatesName,fPfCandidates,fPfCandidatesFromBranch);
  if (!fPfCandidates) {
    SendError(kAbortModule,"Process","Pointer to input PF Cands collection %s null.",fPfCandidatesName.Data());
    return;
  }

  // Create output collections
  fOutputJets = new JetOArr;
  fOutputJets->SetName(fOutputJetsName);
  fOutputJets->SetOwner(kTRUE);
    
  std::vector<fastjet::PseudoJet> fjParts;
  // Push all particle flow candidates of the input PFjet into fastjet particle collection
  for (UInt_t j=0; j<fPfCandidates->GetEntries(); ++j) {
    const PFCandidate *pfCand = fPfCandidates->At(j);
    // Exclude very soft (unphysical) particles
    if (pfCand->Pt() < fParticleMinPt)
      continue;
    fjParts.push_back(fastjet::PseudoJet(pfCand->Px(),pfCand->Py(),pfCand->Pz(),pfCand->E()));
    fjParts.back().set_user_index(j);
  }	
  
  // Setup the clusters for fastjet
  fastjet::ClusterSequenceArea *fjClustering =
    new fastjet::ClusterSequenceArea(fjParts,*fJetDef,*fAreaDefinition);

  // ---- Fastjet is ready ----

  // Produce a new set of jets based on the fastjet particle collection and the defined clustering
  // Cut off fat jets with pt < fJetMinPt GeV
  std::vector<fastjet::PseudoJet> fjOutJets = sorted_by_pt(fjClustering->inclusive_jets(fJetMinPt)); 
  // Check that the output collection size is non-null, otherwise nothing to be done further
  if (fjOutJets.size() < 1) {
    printf(" FastJetMod - WARNING - input PFCands produces null reclustering output!\n");

    if ((fjClustering->inclusive_jets()).size() > 0) 
      fjClustering->delete_self_when_unused();
    delete fjClustering;

    return;
  }

  // Now loop over PFJets and fill the output collection
  for (UInt_t j=0; j<fjOutJets.size(); ++j) {
    // Inizialize PFJet with 4-vector
    PFJet* outJet = new PFJet(fjOutJets[j].px(),
                              fjOutJets[j].py(),
                              fjOutJets[j].pz(),
                              fjOutJets[j].e());

    // Setup PFJet area
    outJet->SetJetArea(fjOutJets[j].area());
    
    // Setup PFJet particle flow related quantities
    FillPFJet(outJet, fjOutJets[j]);
                              
    // Add this jet to the output collection
    fOutputJets->AddOwned(outJet);
        
  } //end loop on fastjet jets
  
  // Now sort the output collections
  fOutputJets->Sort();
      
  // add to event for other modules to use
  AddObjThisEvt(fOutputJets);  
  
  // some memory cleanup
  fjClustering->delete_self_when_unused();
  delete fjClustering;
  
  return;
}

//--------------------------------------------------------------------------------------------------
void FastJetMod::SlaveBegin()
{
  // Run startup code on the computer (slave) doing the actual analysis. 
  if (fGetMatchBtag)
    ReqEventObject(fBtaggedJetsName,fBtaggedJets,fBtaggedJetsFromBranch);
  if (fUseBambuJets)
    ReqEventObject(fJetsName,fJets,fJetsFromBranch);
  ReqEventObject(fPfCandidatesName,fPfCandidates,fPfCandidatesFromBranch);
    
  // AKT constructor
  fJetDef = new fastjet::JetDefinition(fastjet::cambridge_algorithm, fJetConeSize);
  
  // Initialize area caculation (done with ghost particles)
  int activeAreaRepeats = 1;
  double ghostArea = 0.01;
  double ghostEtaMax = 7.0;
  fActiveArea = new fastjet::GhostedAreaSpec(ghostEtaMax,activeAreaRepeats,ghostArea);
  fAreaDefinition = new fastjet::AreaDefinition(fastjet::active_area_explicit_ghosts,*fActiveArea);  
    
  return;
}

//--------------------------------------------------------------------------------------------------
void FastJetMod::SlaveTerminate()
{
}

//--------------------------------------------------------------------------------------------------
void FastJetMod::FillPFJet (PFJet *pPFJet, fastjet::PseudoJet &fjJet)
{
  // Prepare the jet observables related to PFConstituents
  float chargedHadronEnergy = 0.;
  float neutralHadronEnergy = 0.;
  float chargedEmEnergy     = 0.;
  float chargedMuEnergy     = 0.;
  float neutralEmEnergy     = 0.;
  int chargedMultiplicity   = 0;
  int neutralMultiplicity   = 0;
  int muonMultiplicity      = 0;
  
  // Loop on input jet constituents vector and discard very soft particles (ghosts)
  for (unsigned int iPart = 0; iPart < fjJet.constituents().size(); iPart++) {
    if (fjJet.constituents()[iPart].perp() < fParticleMinPt)
      continue;
    int thisPFCandIndex = fjJet.constituents()[iPart].user_index();
    // First of all fix the linking between PFJets and PFCandidates
    const PFCandidate *pfCand = fPfCandidates->At(thisPFCandIndex);
    // Check that the pfCandidate exists
    if (!pfCand) {
      printf(" FastJetMod::FillPFJet - WARNING - input PFCand pointer is null, skipping this candidate.");
      continue;
    }    
    
    pPFJet->AddPFCand(pfCand);      
    
    // Now take care of energy fraction and multiplicities
    if (pfCand->PFType() == PFCandidate::eHadron) {
      chargedHadronEnergy += pfCand->E();
      chargedMultiplicity ++;
    }
    if (pfCand->PFType() == PFCandidate::eNeutralHadron) {
      neutralHadronEnergy += pfCand->E();
      neutralMultiplicity ++;
    }
    if (pfCand->PFType() == PFCandidate::eGamma) {
      neutralEmEnergy += pfCand->E();
      neutralMultiplicity ++;
    }
    if (pfCand->PFType() == PFCandidate::eElectron) {
      chargedEmEnergy += pfCand->E();
      chargedMultiplicity ++;
    }
    if (pfCand->PFType() == PFCandidate::eMuon) {
      chargedMuEnergy += pfCand->E();
      chargedMultiplicity ++;
    }    
  }// end loop on jet constituents
  
  // Fill in the energy fractions and multiplicieties
  pPFJet->SetChargedHadronEnergy(chargedHadronEnergy);
  pPFJet->SetNeutralHadronEnergy(neutralHadronEnergy);
  pPFJet->SetChargedEmEnergy(chargedEmEnergy);
  pPFJet->SetChargedMuEnergy(chargedMuEnergy);
  pPFJet->SetNeutralEmEnergy(neutralEmEnergy);
  pPFJet->SetChargedMultiplicity(chargedMultiplicity);
  pPFJet->SetNeutralMultiplicity(neutralMultiplicity);
  pPFJet->SetMuonMultiplicity(muonMultiplicity);   
  
  return;
}

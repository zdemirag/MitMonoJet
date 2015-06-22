#include "MitMonoJet/Mods/interface/PreRun2SynchExercise.h"

#include "MitAna/DataTree/interface/VertexCol.h"
#include "MitAna/DataTree/interface/PFMetCol.h"
#include "MitAna/DataTree/interface/PFJetCol.h"
#include "MitAna/DataTree/interface/ElectronCol.h"
#include "MitAna/DataTree/interface/MuonCol.h"
#include "MitAna/DataTree/interface/PFTauCol.h"
#include "MitAna/DataTree/interface/PhotonCol.h"

#include "TVector2.h"

ClassImp(mithep::PreRun2SynchExercise)

void
mithep::PreRun2SynchExercise::Process()
{
  auto* inVertices = GetObject<mithep::VertexCol>(fVerticesName);
  auto* inMets = GetObject<mithep::PFMetCol>(fMetName);
  auto* inJets = GetObject<mithep::PFJetCol>(fJetsName);
  auto* inElectrons = GetObject<mithep::ElectronCol>(fElectronsName);
  auto* inMuons = GetObject<mithep::MuonCol>(fMuonsName);
  auto* inTaus = GetObject<mithep::PFTauCol>(fTausName);
  auto* inPhotons = GetObject<mithep::PhotonCol>(fPhotonsName);

  fSynchPass[kGoodVertex] = inVertices->GetEntries() != 0;
  fSynchPass[kMET200] = inMets->At(0)->Pt() > 200.;
  fSynchPass[kNAK4Jets] = inJets->GetEntries() < 3;
  fSynchPass[kEMuVeto] = inElectrons->GetEntries() + inMuons->GetEntries() == 0;
  fSynchPass[kTauVeto] = inTaus->GetEntries() == 0;
  fSynchPass[kPhotonVeto] = inPhotons->GetEntries() == 0;

  PFJet const* leadJet = 0;
  PFJet const* trailJet = 0;

  fSynchPass[kCleanJets] = false;
  fSynchPass[kLeadingJet110] = false;
  fSynchPass[kDeltaPhiJ1J2] = false;
  if (inJets->GetEntries() > 0) {
    leadJet = inJets->At(0);
    if (inJets->GetEntries() == 2) {
      trailJet = inJets->At(1);
      if (inJets->At(1)->Pt() > leadJet->Pt()) {
        leadJet = inJets->At(1);
        trailJet = inJets->At(0);
      }
    }

    if (leadJet->ChargedHadronEnergy() / leadJet->E() > 0.2 &&
        leadJet->NeutralHadronEnergy() / leadJet->E() < 0.7 &&
        leadJet->NeutralEmEnergy() / leadJet->E() < 0.7) {
      fSynchPass[kCleanJets] = true;
      if (leadJet->Pt() > 110.)
        fSynchPass[kLeadingJet110] = true;
    }

    fSynchPass[kDeltaPhiJ1J2] = true;

    if (trailJet) {
      if (trailJet->NeutralHadronEnergy() / trailJet->E() < 0.7 &&
          trailJet->NeutralEmEnergy() / trailJet->E() < 0.9) {
        if (std::abs(TVector2::Phi_mpi_pi(leadJet->Phi() - trailJet->Phi())) > 2.5)
          fSynchPass[kDeltaPhiJ1J2] = false;
      }
      else {
        fSynchPass[kCleanJets] = false;
        fSynchPass[kDeltaPhiJ1J2] = false;
      }
    }
  }

  fTree->Fill();
}

void
mithep::PreRun2SynchExercise::SlaveBegin()
{
  fTree = new TTree("synchTree", "Synch checkpoints");
  fTree->Branch("GoodVertex", fSynchPass + kGoodVertex, "pass/O");
  fTree->Branch("LeadingJet110", fSynchPass + kLeadingJet110, "pass/O");
  fTree->Branch("DeltaPhiJ1J2", fSynchPass + kDeltaPhiJ1J2, "pass/O"); 
  fTree->Branch("MET200", fSynchPass + kMET200, "pass/O");       
  fTree->Branch("NAK4Jets", fSynchPass + kNAK4Jets, "pass/O");     
  fTree->Branch("EMuVeto", fSynchPass + kEMuVeto, "pass/O");      
  fTree->Branch("TauVeto", fSynchPass + kTauVeto, "pass/O");      
  fTree->Branch("PhotonVeto", fSynchPass + kPhotonVeto, "pass/O");   

  AddOutput(fTree);
}

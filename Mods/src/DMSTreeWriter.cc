#include <TSystem.h>
#include <TFile.h>
#include "MitCommon/MathTools/interface/MathUtils.h"
#include "MitAna/DataTree/interface/MuonFwd.h"
#include "MitAna/DataTree/interface/ElectronFwd.h"
#include "MitAna/DataTree/interface/VertexCol.h"
#include "MitAna/TreeMod/interface/BaseMod.h"
#include "MitAna/DataTree/interface/PhotonCol.h"
#include "MitAna/DataTree/interface/GenericParticle.h"
#include "MitAna/DataTree/interface/PFMet.h"
#include "MitAna/DataTree/interface/PFJetCol.h"
#include "MitAna/DataTree/interface/PFTauCol.h"
#include "MitPhysics/Init/interface/ModNames.h"
#include "MitPhysics/Utils/interface/MuonTools.h"
#include "MitMonoJet/DataTree/interface/XlFatJet.h"

#include "MitMonoJet/Mods/interface/DMSTreeWriter.h"

using namespace mithep;

ClassImp(mithep::DMSTreeWriter)

//--------------------------------------------------------------------------------------------------
DMSTreeWriter::DMSTreeWriter(const char *name, const char *title) :
  BaseMod                 (name,title),
  fEvtSelDataName         (Names::gkEvtSelDataBrn),
  fRawMetName             ("PFMet"),
  fMetName                ("PFMet"),
  fMetMVAName             ("PFMetMVA"),
  fPhotonsName            (Names::gkPhotonBrn),
  fElectronsName          (Names::gkElectronBrn),
  fMuonsName              (Names::gkMuonBrn),
  fTausName               (Names::gkPFTauBrn),
  fJetsName               (Names::gkPFJetBrn),
  fFatJetsName            ("XlFatJets"),
  fSubJetsName            ("XlSubJets"),
  fPVName                 (Names::gkPVBeamSpotBrn),
  fPileUpDenName          (Names::gkPileupEnergyDensityBrn),
  fPileUpName             (Names::gkPileupInfoBrn),
  fTriggerObjectsName     ("MyHltPhotObjs"),
  fIsData                 (false),
  fMetFromBranch          (kTRUE),
  fMetMVAFromBranch       (kTRUE),
  fPhotonsFromBranch      (kTRUE),
  fElectronsFromBranch    (kTRUE),
  fMuonsFromBranch        (kTRUE),
  fTausFromBranch         (kTRUE),
  fJetsFromBranch         (kTRUE),
  fFatJetsFromBranch      (kTRUE),
  fSubJetsFromBranch      (kTRUE),
  fPVFromBranch           (kTRUE),
  // -------------------------
  fRawMet                 (0),
  fMet                    (0),
  fMetMVA                 (0),
  fPhotons                (0),
  fElectrons              (0),
  fMuons                  (0),
  fPFTaus                 (0),
  fJets                   (0),
  fFatJets                (0),
  fSubJets                (0),
  fPV                     (0),
  fPileUp                 (0),
  fPileUpDen              (0),
  fEvtSelData             (0),
  fTrigObj                (0),
  // -------------------------
  fOutputFile             (0)

{
  // Constructor
}

DMSTreeWriter::~DMSTreeWriter()
{
  // Destructor
  fOutputFile->Close();
}

//--------------------------------------------------------------------------------------------------
void DMSTreeWriter::SlaveTerminate()
{
  fOutputFile->WriteTObject(fMitDMSTree.tree_,fMitDMSTree.tree_->GetName());
}

//--------------------------------------------------------------------------------------------------
void DMSTreeWriter::Process()
{
  // Process entries of the tree.
  LoadEventObject(fEvtSelDataName,    fEvtSelData,    true);
  LoadEventObject(fPileUpDenName,     fPileUpDen,     true);
  if (!fIsData) {
    LoadBranch(fPileUpName);
  }
  LoadEventObject(fTriggerObjectsName,fTrigObj,       true);

  LoadEventObject(fRawMetName,        fRawMet,        true);
  LoadEventObject(fMetName,           fMet,           fMetFromBranch);
  LoadEventObject(fMetMVAName,        fMetMVA,        fMetMVAFromBranch);
  LoadEventObject(fPhotonsName,       fPhotons,       fPhotonsFromBranch);
  LoadEventObject(fElectronsName,     fElectrons,     fElectronsFromBranch);
  LoadEventObject(fMuonsName,         fMuons,         fMuonsFromBranch);
  LoadEventObject(fTausName,          fPFTaus,        fTausFromBranch);
  LoadEventObject(fJetsName,          fJets,          fJetsFromBranch);
  LoadEventObject(fFatJetsName,       fFatJets,       fJetsFromBranch);
  LoadEventObject(fSubJetsName,       fSubJets,       fJetsFromBranch);
  LoadEventObject(fPVName,            fPV,            fPVFromBranch);

  // initialize the tree variables
  fMitDMSTree.InitVariables();

  // EVTSELDATA
  fMitDMSTree.metFiltersWord_ = fEvtSelData->metFiltersWord();

  // PILEUP RELATED

  if (! fIsData) {
    // loop over the pileup summary info and grab what you need
    for (UInt_t i=0; i<fPileUp->GetEntries(); ++i) {
      const PileupInfo *puinfo = fPileUp->At(i);
      if (puinfo->GetBunchCrossing() ==  0)
        fMitDMSTree.npu_ = puinfo->GetPU_NumMean();
      if (puinfo->GetBunchCrossing() ==  1)
        fMitDMSTree.npuPlusOne_ = puinfo->GetPU_NumInteractions();
      if (puinfo->GetBunchCrossing() == -1)
        fMitDMSTree.npuMinusOne_ = puinfo->GetPU_NumInteractions();
    }
  }

  // EVENT

  fMitDMSTree.run_   = GetEventHeader()->RunNum();
  fMitDMSTree.lumi_  = GetEventHeader()->LumiSec();
  fMitDMSTree.event_ = GetEventHeader()->EvtNum();
  fMitDMSTree.nvtx_  = fPV->GetEntries();

  // HLT
  fMitDMSTree.trigger_ = 0;
  
  if (! fTrigObj)
    printf("MonoJetTreeWriter::TriggerObjectCol not found\n");

  else {
    // loop through the stored trigger objects and find corresponding trigger name
    int nGoodCntJets = 0;
    bool hasGoodMET = 0;
    bool hasGoodMHT = 0;
    bool hasGoodMuons = 0;
    for (UInt_t i=0;i<fTrigObj->GetEntries();++i) {
      const TriggerObject *to = fTrigObj->At(i);
      if (to->TriggerType() == 83)
        hasGoodMuons = true;
      if (to->TriggerType() == 85 && to->Pt() > 80 && fabs(to->Eta()) < 2.4)
        nGoodCntJets++;
      if (to->TriggerType() == 87)
        hasGoodMET = true;
      if (to->TriggerType() == 90)
        hasGoodMHT = true;
    }
    // default MonoJet
    if (nGoodCntJets > 0 && hasGoodMHT)
      fMitDMSTree.trigger_ |= 1 << 0;
    if (hasGoodMET)
      fMitDMSTree.trigger_ |= 1 << 1;
    // default single muon
    if (hasGoodMuons)
      fMitDMSTree.trigger_ |= 1 << 2;
  }

  // MET BASICS

  fMitDMSTree.metRaw_        = fRawMet->At(0)->Pt();
  fMitDMSTree.metRawPhi_     = fRawMet->At(0)->Phi();
  fMitDMSTree.met_           = fMet->At(0)->Pt();
  fMitDMSTree.metPhi_        = fMet->At(0)->Phi();
  fMitDMSTree.mvamet_        = fMetMVA->At(0)->Pt();
  fMitDMSTree.mvametPhi_     = fMetMVA->At(0)->Phi();
  fMitDMSTree.mvaCov00_      = fMetMVA->At(0)->Cov00();
  fMitDMSTree.mvaCov10_      = fMetMVA->At(0)->Cov10();
  fMitDMSTree.mvaCov01_      = fMetMVA->At(0)->Cov01();
  fMitDMSTree.mvaCov11_      = fMetMVA->At(0)->Cov11();

  // LEPTONS (MU+ELE)

  fMitDMSTree.nlep_ = fMuons->GetEntries() + fElectrons->GetEntries();
  if (fMuons->GetEntries() > 1) {
    fMitDMSTree.lep1_ = fMuons->At(0)->Mom();
    fMitDMSTree.lid1_ = 13;
    fMitDMSTree.lep2_ = fMuons->At(1)->Mom();
    fMitDMSTree.lid2_ = 13;
  }
  else if (fMuons->GetEntries() > 0) {
    fMitDMSTree.lep1_ = fMuons->At(0)->Mom();
    fMitDMSTree.lid1_ = 13;
    //mu e
    if (fElectrons->GetEntries() > 0) {
      fMitDMSTree.lep2_ = fElectrons->At(0)->Mom();
      fMitDMSTree.lid2_ = 11;
    }             
  }
  else if (fElectrons->GetEntries() > 1) {
    //e e
    fMitDMSTree.lep1_ = fElectrons->At(0)->Mom();
    fMitDMSTree.lid1_ = 13;
    fMitDMSTree.lep2_ = fElectrons->At(1)->Mom();
    fMitDMSTree.lid2_ = 13;
  }
  else if (fElectrons->GetEntries() > 0) {
    //e e
    fMitDMSTree.lep1_ = fElectrons->At(0)->Mom();
    fMitDMSTree.lid1_ = 13;
  }    

  // TAUS

  fMitDMSTree.ntaus_ = fPFTaus->GetEntries();
  if (fPFTaus->GetEntries() >= 1) {
    fMitDMSTree.tau1_ = fPFTaus->At(0)->Mom();
  }

  // PHOTON(S)

  fMitDMSTree.nphotons_ = fPhotons->GetEntries();
  if (fPhotons->GetEntries() >= 1) {
    const Photon *photon = fPhotons->At(0);
    fMitDMSTree.pho1_ = photon->Mom();
  }

  // FAT JET (select highest in pt)
  if (fFatJets->GetEntries() >= 1) {
    const XlFatJet *fjet = fFatJets->At(0);    
    fMitDMSTree.tjet_        = fjet->Mom();
    fMitDMSTree.tjetGroomed_ = fjet->GroomedMom();
    fMitDMSTree.tjetBtag_    = fjet->CombinedSecondaryVertexBJetTagsDisc();
    fMitDMSTree.tjetTau1_    = fjet->Tau1();
    fMitDMSTree.tjetTau2_    = fjet->Tau2();
    fMitDMSTree.tjetTau3_    = fjet->Tau3();
    fMitDMSTree.tjetC2b0_    = fjet->C2b0();
    fMitDMSTree.tjetPartonId_= -1;
  
    fMitDMSTree.nsjets_ = fjet->NSubJets();
    if (fMitDMSTree.nsjets_ >= 1) {
      fMitDMSTree.sjet1_        = fjet->SubJet(0)->Mom();
      fMitDMSTree.sjet2_        = fjet->SubJet(1)->Mom();
    }
  }
  
  // JETS : careful since the hardest overlaps with the tag jet
  fMitDMSTree.njets_ = fJets->GetEntries();

  // PILEUP RELATED

  if (! fIsData) {
    // loop over the pileup summary info and grab what you need
    for (UInt_t i=0; i<fPileUp->GetEntries(); ++i) {
      const PileupInfo *puinfo = fPileUp->At(i);
      if (puinfo->GetBunchCrossing() ==  0)
        fMitDMSTree.npu_ = puinfo->GetPU_NumMean();
      if (puinfo->GetBunchCrossing() ==  1)
        fMitDMSTree.npuPlusOne_ = puinfo->GetPU_NumInteractions();
      if (puinfo->GetBunchCrossing() == -1)
        fMitDMSTree.npuMinusOne_ = puinfo->GetPU_NumInteractions();
    }
  }

  // Finally fill the tree
  fMitDMSTree.tree_->Fill();

  return;
}

//--------------------------------------------------------------------------------------------------
void DMSTreeWriter::SlaveBegin()
{
  // Run startup code on the computer (slave) doing the actual analysis. Here, we request all
  // relevant information.

  if (! fIsData) {
    ReqBranch(fPileUpName,           fPileUp);
  }
  ReqEventObject(fPileUpDenName,     fPileUpDen,     true);
  ReqEventObject(fPVName,            fPV,            fPVFromBranch);
  ReqEventObject("EvtSelData",       fEvtSelData,    true);
  ReqEventObject(fTriggerObjectsName,fTrigObj,       true);

  ReqEventObject(fPhotonsName,       fPhotons,       fPhotonsFromBranch);
  ReqEventObject(fElectronsName,     fElectrons,     fElectronsFromBranch);
  ReqEventObject(fMuonsName,         fMuons,         fMuonsFromBranch);
  ReqEventObject(fTausName,          fPFTaus,        fTausFromBranch);
  ReqEventObject(fJetsName,          fJets,          fJetsFromBranch);
  ReqEventObject(fFatJetsName,       fFatJets,       fFatJetsFromBranch);
  ReqEventObject(fSubJetsName,       fSubJets,       fSubJetsFromBranch);
  ReqEventObject(fRawMetName,        fRawMet,        true);
  ReqEventObject(fMetName,           fMet,           fMetFromBranch);
  ReqEventObject(fMetMVAName,        fMetMVA,        fMetMVAFromBranch);

  // Create Ntuple Tree
  fOutputFile = TFile::Open(TString::Format("%s_tmp.root",GetName()),"RECREATE");
  fMitDMSTree.CreateTree(0);
  fMitDMSTree.tree_->SetAutoSave(300e9);
  fMitDMSTree.tree_->SetDirectory(fOutputFile);
  AddOutput(fMitDMSTree.tree_);
}

void DMSTreeWriter::CorrectMet(const float met, const float metPhi,
				   const Particle *l1, const Particle *l2,
                                   float &newMet, float &newMetPhi)
{
  // inputs:  met, metPhi, l1,  [ l2  only used if pointer is non-zero ]
  // outputs: newMet, newMetPhi

  float newMetX;
  float newMetY;
  if (l2) {   // these are doubles on the right: need full calculation in one line
    newMetX = met*TMath::Cos(metPhi) + l1->Mom().Px() + l2->Mom().Px();
    newMetY = met*TMath::Sin(metPhi) + l1->Mom().Py() + l2->Mom().Py();
  }
  else {
    newMetX = met*TMath::Cos(metPhi) + l1->Mom().Px();
    newMetY = met*TMath::Sin(metPhi) + l1->Mom().Py();
  }
  
  newMet    = TMath::Sqrt(TMath::Power(newMetX,2) + TMath::Power(newMetY,2));
  newMetPhi = TMath::ATan2(newMetY,newMetX);
}

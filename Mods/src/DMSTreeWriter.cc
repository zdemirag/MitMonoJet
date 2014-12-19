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
#include "MitMonoJet/DataTree/interface/XsIsoParticle.h"

#include "MitMonoJet/Mods/interface/DMSTreeWriter.h"

using namespace mithep;

ClassImp(mithep::DMSTreeWriter)

//--------------------------------------------------------------------------------------------------
DMSTreeWriter::DMSTreeWriter(const char *name, const char *title) :
  BaseMod                 (name,title),
  fEvtSelDataName         ("XlEvtSelData"),
  fRawMetName             ("PFMet"),
  fMetMVAName             ("PFMetMVA"),
  fPhotonsName            ("XsPhotons"),
  fElectronsName          ("XsElectrons"),
  fMuonsName              ("XsMuons"),
  fTausName               ("XsTaus"),
  fJetsName               (Names::gkPFJetBrn),
  fFatJetsName            ("XlFatJets"),
  fSubJetsName            ("XlSubJets"),
  fPVName                 (Names::gkPVBeamSpotBrn),
  fPileUpDenName          (Names::gkPileupEnergyDensityBrn),
  fPileUpName             (Names::gkPileupInfoBrn),
  fMCEventInfoName        (Names::gkMCEvtInfoBrn),
  fMCParticlesName        (Names::gkMCPartBrn),
  fTriggerObjectsName     ("MyHltPhotObjs"),
  fIsData                 (false),
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
  fMCEventInfo            (0),
  fMCParticles            (0),
  fEvtSelData             (0),
  fTrigObj                (0),
  fPUInputFileName        ("MyInputPUFile"),       
  fPUTargetFileName       ("MyTargetPUFile"),
  fPUInput                (0),       
  fPUTarget               (0),
  fPUWeight               (0),
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
    LoadBranch(fMCEventInfoName);
    LoadBranch(fMCParticlesName);
  }
  LoadEventObject(fTriggerObjectsName,fTrigObj,       true);

  LoadEventObject(fRawMetName,        fRawMet,        true);
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
  fMitDMSTree.preselWord_     = fEvtSelData->preselWord();

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
    bool hasGoodPhotons = 0;
    for (UInt_t i=0;i<fTrigObj->GetEntries();++i) {
      const TriggerObject *to = fTrigObj->At(i);
      //to->Print(); 
      if (to->TriggerType() == TriggerObject::TriggerJet 
       && to->Pt() > 80 && fabs(to->Eta()) < 2.4)
        nGoodCntJets++;
      if (to->TriggerType() == TriggerObject::TriggerMHT)
        hasGoodMHT = true;
      if (to->TriggerType() == TriggerObject::TriggerMET)
        hasGoodMET = true;
      if (to->TriggerType() == TriggerObject::TriggerMuon 
       && to->Pt() > 24 && fabs(to->Eta()) < 2.1)
        hasGoodMuons = true;
      if (to->TriggerType() == TriggerObject::TriggerPhoton
       && to->Pt() > 130)
        hasGoodPhotons = true;
    }
    // default MonoJet
    if (nGoodCntJets > 0 && hasGoodMHT)
      fMitDMSTree.trigger_ |= MitDMSTree::HLTJetMet;
    if (hasGoodMET)
      fMitDMSTree.trigger_ |= MitDMSTree::HLTMet;
    // default single muon
    if (hasGoodMuons)
      fMitDMSTree.trigger_ |= MitDMSTree::HLTMuon;
    // default single photon
    if (hasGoodPhotons)
      fMitDMSTree.trigger_ |= MitDMSTree::HLTPhoton;
  }

  // MET BASICS

  fMitDMSTree.metRaw_        = fRawMet->At(0)->Pt();
  fMitDMSTree.metRawPhi_     = fRawMet->At(0)->Phi();
  fMitDMSTree.met_           = fMitDMSTree.metRaw_;
  fMitDMSTree.metPhi_        = fMitDMSTree.metRawPhi_;
  fMitDMSTree.mvamet_        = fMetMVA->At(0)->Pt();
  fMitDMSTree.mvametPhi_     = fMetMVA->At(0)->Phi();
  fMitDMSTree.mvaCov00_      = fMetMVA->At(0)->Cov00();
  fMitDMSTree.mvaCov10_      = fMetMVA->At(0)->Cov10();
  fMitDMSTree.mvaCov01_      = fMetMVA->At(0)->Cov01();
  fMitDMSTree.mvaCov11_      = fMetMVA->At(0)->Cov11();

  // LEPTONS (MU+ELE), save tight id for further studies
  // Also perform met/mt/mll computation according to the relevant case (muons only)
  fMitDMSTree.nlep_ = fMuons->GetEntries() + fElectrons->GetEntries();
  if (fMuons->GetEntries() > 1) {
    //mu mu
    fMitDMSTree.lep1_ = fMuons->At(0)->Mom();
    fMitDMSTree.lid1_ = fMuons->At(0)->Charge()*13;
    if (fMuons->At(0)->ParticleId() == XsIsoParticle::EParticleId::eTightMuon)
      fMitDMSTree.lid1_ = fMitDMSTree.lid1_*10;
    if (fMuons->At(0)->ParticleId() == XsIsoParticle::EParticleId::eIsoMuon)
      fMitDMSTree.lid1_ = fMitDMSTree.lid1_*100;

    fMitDMSTree.lep2_ = fMuons->At(1)->Mom();
    fMitDMSTree.lid2_ = fMuons->At(1)->Charge()*13;
    if (fMuons->At(1)->ParticleId() == XsIsoParticle::EParticleId::eTightMuon)
      fMitDMSTree.lid2_ = fMitDMSTree.lid2_*10;
    if (fMuons->At(1)->ParticleId() == XsIsoParticle::EParticleId::eIsoMuon)
      fMitDMSTree.lid2_ = fMitDMSTree.lid2_*100;

    CorrectMet(fMitDMSTree.metRaw_,fMitDMSTree.metRawPhi_,fMitDMSTree.lep1_,fMitDMSTree.lep2_,
               fMitDMSTree.met_,fMitDMSTree.metPhi_);    
    fMitDMSTree.mll_ = (fMitDMSTree.lep1_ + fMitDMSTree.lep2_).M();
  }
  else if (fMuons->GetEntries() == 1) {
    fMitDMSTree.lep1_ = fMuons->At(0)->Mom();
    fMitDMSTree.lid1_ = fMuons->At(0)->Charge()*13;
    if (fMuons->At(0)->ParticleId() == XsIsoParticle::EParticleId::eTightMuon)
      fMitDMSTree.lid1_ = fMitDMSTree.lid1_*10;
    if (fMuons->At(0)->ParticleId() == XsIsoParticle::EParticleId::eIsoMuon)
      fMitDMSTree.lid1_ = fMitDMSTree.lid1_*100;

    CorrectMet(fMitDMSTree.metRaw_,fMitDMSTree.metRawPhi_,fMitDMSTree.lep1_,
               fMitDMSTree.met_,fMitDMSTree.metPhi_);    
    fMitDMSTree.mt_ = GetMt(fMitDMSTree.lep1_,fMitDMSTree.metRaw_,fMitDMSTree.metRawPhi_);
    //mu e
    if (fElectrons->GetEntries() > 0) {
      fMitDMSTree.lep2_ = fElectrons->At(0)->Mom();
      fMitDMSTree.lid2_ = fElectrons->At(0)->Charge()*11;
    }             
  }
  //e e
  else if (fMuons->GetEntries() == 0 && fElectrons->GetEntries() > 1) {
    fMitDMSTree.lep1_ = fElectrons->At(0)->Mom();
    fMitDMSTree.lid1_ = fElectrons->At(0)->Charge()*11;
    fMitDMSTree.lep2_ = fElectrons->At(1)->Mom();
    fMitDMSTree.lid2_ = fElectrons->At(1)->Charge()*11;
  }
  else if (fMuons->GetEntries() == 0 && fElectrons->GetEntries() > 0) {
    //e
    fMitDMSTree.lep1_ = fElectrons->At(0)->Mom();
    fMitDMSTree.lid1_ = fElectrons->At(0)->Charge()*11;
  }    
  // Make Muon-HLT matching 
  if (fMitDMSTree.nlep_ && fMitDMSTree.lid1_ == 13
   && IsHLTMatched(fMitDMSTree.lep1_, 
                   TriggerObject::TriggerMuon,
                   0.3, 
                   24., 2.1))
    fMitDMSTree.HLTmatch_ |= MitDMSTree::MuonMatch;         

  // TAUS

  fMitDMSTree.ntaus_ = fPFTaus->GetEntries();
  if (fPFTaus->GetEntries() >= 1) {
    fMitDMSTree.tau1_ = fPFTaus->At(0)->Mom();
  }

  // PHOTON(S)

  fMitDMSTree.nphotons_ = fPhotons->GetEntries();
  if (fPhotons->GetEntries() >= 1) {
    const XsIsoParticle *photon = fPhotons->At(0);
    fMitDMSTree.pho1_ = photon->Mom();
    // Make Photon-HLT matching 
    if (IsHLTMatched(fMitDMSTree.pho1_, 
                     TriggerObject::TriggerPhoton,
                     0.3, 
                     130.))
      fMitDMSTree.HLTmatch_ |= MitDMSTree::PhotonMatch;         
    // Correct MET (only if no leptons and high pt photon!)
    if (fMitDMSTree.nlep_ == 0 && photon->Pt() > 150.) 
      CorrectMet(fMitDMSTree.metRaw_,fMitDMSTree.metRawPhi_,fMitDMSTree.pho1_,
                 fMitDMSTree.met_,fMitDMSTree.metPhi_);    
  }

  // FAT JETS  
  fMitDMSTree.nfjets_ = fFatJets->GetEntries();
  for (UInt_t i = 0; i < fFatJets->GetEntries(); ++i) {
    
    if (i == 0) {
      const XlFatJet *fjet = fFatJets->At(i);    
      fMitDMSTree.fjet1_       = fjet->Mom();
      // Further cleaning for larger cones
      if (!fjetIsCleaned(fMitDMSTree.fjet1_,0.5))
        continue;

      fMitDMSTree.fjet1Btag_    = GetFatJetBtag(fMitDMSTree.fjet1_, 0.5);
      fMitDMSTree.fjet1Charge_  = fjet->Charge();
      fMitDMSTree.fjet1QGtag_   = fjet->QGTag();
      fMitDMSTree.fjet1Tau1_    = fjet->Tau1();
      fMitDMSTree.fjet1Tau2_    = fjet->Tau2();
      fMitDMSTree.fjet1Tau3_    = fjet->Tau3();
      fMitDMSTree.fjet1C2b0_    = fjet->C2b0();
      fMitDMSTree.fjet1C2b0p2_      = fjet->C2b0p2();      
      fMitDMSTree.fjet1C2b0p5_      = fjet->C2b0p5();      
      fMitDMSTree.fjet1C2b1_        = fjet->C2b1();        
      fMitDMSTree.fjet1C2b2_        = fjet->C2b2();        
      fMitDMSTree.fjet1QJetVol_     = fjet->QJetVol();     
      fMitDMSTree.fjet1MassSDbm1_   = fjet->MassSDbm1();   
      fMitDMSTree.fjet1MassSDb0_    = fjet->MassSDb0();    
      fMitDMSTree.fjet1MassSDb1_    = fjet->MassSDb1();    
      fMitDMSTree.fjet1MassSDb2_    = fjet->MassSDb2();    
      fMitDMSTree.fjet1MassPruned_  = fjet->MassPruned();  
      fMitDMSTree.fjet1MassFiltered_= fjet->MassFiltered();
      fMitDMSTree.fjet1MassTrimmed_ = fjet->MassTrimmed();
      fMitDMSTree.fjet1Pull_        = fjet->Pull();
      fMitDMSTree.fjet1PullAngle_   = fjet->PullAngle();
      if (!fIsData)  
        fMitDMSTree.fjet1PartonId_  = JetPartonMatch(fMitDMSTree.fjet1_, 0.7);  
    
      fMitDMSTree.fjet1nsj_ = fjet->NSubJets();
      if (fMitDMSTree.fjet1nsj_ > 0) {
        fMitDMSTree.fjet1sj1_         = fjet->SubJet(0)->Mom();
        fMitDMSTree.fjet1QGtagSub1_   = fjet->SubJet(0)->QGTag();
        fMitDMSTree.fjet1QGPtDSub1_   = fjet->SubJet(0)->QGPtD();
        fMitDMSTree.fjet1QGAxis1Sub1_ = fjet->SubJet(0)->QGAxis1();
        fMitDMSTree.fjet1QGAxis2Sub1_ = fjet->SubJet(0)->QGAxis2();
        fMitDMSTree.fjet1QGMultSub1_  = fjet->SubJet(0)->QGMult();
      }
      if (fMitDMSTree.fjet1nsj_ > 1) {
        fMitDMSTree.fjet1sj2_         = fjet->SubJet(1)->Mom();
        fMitDMSTree.fjet1QGtagSub2_   = fjet->SubJet(1)->QGTag();
        fMitDMSTree.fjet1QGPtDSub2_   = fjet->SubJet(1)->QGPtD();
        fMitDMSTree.fjet1QGAxis1Sub2_ = fjet->SubJet(1)->QGAxis1();
        fMitDMSTree.fjet1QGAxis2Sub2_ = fjet->SubJet(1)->QGAxis2();
        fMitDMSTree.fjet1QGMultSub2_  = fjet->SubJet(1)->QGMult();
      }

    }// end filling of first fat jet

    if (i == 1) {
      const XlFatJet *fjet = fFatJets->At(i);    
      fMitDMSTree.fjet2_       = fjet->Mom();
      // Further cleaning for larger cones
      if (!fjetIsCleaned(fMitDMSTree.fjet2_,0.5))
        continue;

      fMitDMSTree.fjet2Btag_    = GetFatJetBtag(fMitDMSTree.fjet2_, 0.5);
      fMitDMSTree.fjet2Charge_  = fjet->Charge();
      fMitDMSTree.fjet2QGtag_   = fjet->QGTag();
      fMitDMSTree.fjet2Tau1_    = fjet->Tau1();
      fMitDMSTree.fjet2Tau2_    = fjet->Tau2();
      fMitDMSTree.fjet2Tau3_    = fjet->Tau3();
      fMitDMSTree.fjet2C2b0_    = fjet->C2b0();
      fMitDMSTree.fjet2C2b0p2_      = fjet->C2b0p2();      
      fMitDMSTree.fjet2C2b0p5_      = fjet->C2b0p5();      
      fMitDMSTree.fjet2C2b1_        = fjet->C2b1();        
      fMitDMSTree.fjet2C2b2_        = fjet->C2b2();        
      fMitDMSTree.fjet2QJetVol_     = fjet->QJetVol();     
      fMitDMSTree.fjet2MassSDbm1_   = fjet->MassSDbm1();   
      fMitDMSTree.fjet2MassSDb0_    = fjet->MassSDb0();    
      fMitDMSTree.fjet2MassSDb1_    = fjet->MassSDb1();    
      fMitDMSTree.fjet2MassSDb2_    = fjet->MassSDb2();    
      fMitDMSTree.fjet2MassPruned_  = fjet->MassPruned();  
      fMitDMSTree.fjet2MassFiltered_= fjet->MassFiltered();
      fMitDMSTree.fjet2MassTrimmed_ = fjet->MassTrimmed();
      fMitDMSTree.fjet2Pull_        = fjet->Pull();
      fMitDMSTree.fjet2PullAngle_   = fjet->PullAngle();
      if (!fIsData)  
        fMitDMSTree.fjet2PartonId_  = JetPartonMatch(fMitDMSTree.fjet2_, 0.7);
                
      fMitDMSTree.fjet2nsj_ = fjet->NSubJets();
      
      if (fMitDMSTree.fjet2nsj_ > 0) {
        fMitDMSTree.fjet2sj1_         = fjet->SubJet(0)->Mom();
        fMitDMSTree.fjet2QGtagSub1_   = fjet->SubJet(0)->QGTag();
        fMitDMSTree.fjet2QGPtDSub1_   = fjet->SubJet(0)->QGPtD();
        fMitDMSTree.fjet2QGAxis1Sub1_ = fjet->SubJet(0)->QGAxis1();
        fMitDMSTree.fjet2QGAxis2Sub1_ = fjet->SubJet(0)->QGAxis2();
        fMitDMSTree.fjet2QGMultSub1_  = fjet->SubJet(0)->QGMult();
      }
      if (fMitDMSTree.fjet2nsj_ > 1) {
        fMitDMSTree.fjet2sj2_         = fjet->SubJet(1)->Mom();
        fMitDMSTree.fjet2QGtagSub2_   = fjet->SubJet(1)->QGTag();
        fMitDMSTree.fjet2QGPtDSub2_   = fjet->SubJet(1)->QGPtD();
        fMitDMSTree.fjet2QGAxis1Sub2_ = fjet->SubJet(1)->QGAxis1();
        fMitDMSTree.fjet2QGAxis2Sub2_ = fjet->SubJet(1)->QGAxis2();
        fMitDMSTree.fjet2QGMultSub2_  = fjet->SubJet(1)->QGMult();
      }
        
    }// end filling of second fat jet

  }
  
  // JETS : careful since the hardest could overlap with the fat jets
  fMitDMSTree.njets_ = fJets->GetEntries();
  for (UInt_t i = 0; i < fJets->GetEntries(); ++i) {
    const PFJet *jet = fJets->At(i);

    if (i == 0) {
      fMitDMSTree.jet1_        = jet->Mom();
      fMitDMSTree.jet1CHF_     = jet->ChargedHadronEnergy()/jet->RawMom().E();
      fMitDMSTree.jet1NHF_     = jet->NeutralHadronEnergy()/jet->RawMom().E();
      fMitDMSTree.jet1NEMF_    = jet->NeutralEmEnergy()/jet->RawMom().E();
      // Make Jet-HLT matching: this is used for preselection
      if (IsHLTMatched(fMitDMSTree.jet1_, 
                       TriggerObject::TriggerJet,
                       0.5, 
                       80., 2.4))
        fMitDMSTree.HLTmatch_ |= MitDMSTree::JetMatch;         
    }
    if (i == 1)
      fMitDMSTree.jet2_        = jet->Mom();
    if (i == 2)
      fMitDMSTree.jet3_        = jet->Mom();
    if (i == 3)
      fMitDMSTree.jet4_        = jet->Mom();
    if (i == 4)
      fMitDMSTree.jet5_        = jet->Mom();
  }

  // B-JETS : careful since the hardest could overlap with the fat jets
  fMitDMSTree.nbjets_ = 0;
  for (UInt_t i = 0; i < fJets->GetEntries(); ++i) {
    const Jet *jet = fJets->At(i);
    
    // Check that the jet is b-tagged
    float btag = jet->CombinedSecondaryVertexBJetTagsDisc();  
    if (btag < 0.244)
      continue;

    // Fill the information for the two hardest b-jets
    if (fMitDMSTree.nbjets_ == 0) {
      fMitDMSTree.bjet1_        = jet->Mom();
      fMitDMSTree.bjet1Btag_    = btag;
    }
    if (fMitDMSTree.nbjets_ == 1) {
      fMitDMSTree.bjet2_        = jet->Mom();
      fMitDMSTree.bjet2Btag_    = btag;
    }  

    // Increment the b-jet counter
    fMitDMSTree.nbjets_ ++;

  }

 
  // MC INFORMATION

  Double_t Q = 0.0;
  Int_t    id1 = 0;
  Double_t x1 = 0.0;
  Double_t pdf1 = 0.0;
  Int_t    id2 = 0;
  Double_t x2 = 0.0;
  Double_t pdf2 = 0.0;
  Int_t    processId = 0;
  if (! fIsData) {
    getGenLevelInfo(fMitDMSTree);

    Q         = fMCEventInfo->Scale();
    id1       = fMCEventInfo->Id1();
    x1        = fMCEventInfo->X1();
    pdf1      = fMCEventInfo->Pdf1();
    id2       = fMCEventInfo->Id2();
    x2        = fMCEventInfo->X2();
    pdf2      = fMCEventInfo->Pdf2();
    processId = fMCEventInfo->ProcessId();
  }

  fMitDMSTree.Q_ = Q;
  fMitDMSTree.id1_ = id1;
  fMitDMSTree.x1_ = x1;
  fMitDMSTree.pdf1_ = pdf1;
  fMitDMSTree.id2_ = id2;
  fMitDMSTree.x2_ = x2;
  fMitDMSTree.pdf2_ = pdf2;
  fMitDMSTree.processId_ = processId;

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
    fMitDMSTree.puweight_ = PUWeight(fMitDMSTree.npu_);
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
    ReqBranch(fMCEventInfoName,      fMCEventInfo);
    ReqEventObject(fMCParticlesName, fMCParticles,   true);
  }
  ReqEventObject(fPileUpDenName,     fPileUpDen,     true);
  ReqEventObject(fPVName,            fPV,            fPVFromBranch);
  ReqEventObject(fEvtSelDataName,    fEvtSelData,    true);
  ReqEventObject(fTriggerObjectsName,fTrigObj,       true);

  ReqEventObject(fPhotonsName,       fPhotons,       fPhotonsFromBranch);
  ReqEventObject(fElectronsName,     fElectrons,     fElectronsFromBranch);
  ReqEventObject(fMuonsName,         fMuons,         fMuonsFromBranch);
  ReqEventObject(fTausName,          fPFTaus,        fTausFromBranch);
  ReqEventObject(fJetsName,          fJets,          fJetsFromBranch);
  ReqEventObject(fFatJetsName,       fFatJets,       fFatJetsFromBranch);
  ReqEventObject(fSubJetsName,       fSubJets,       fSubJetsFromBranch);
  ReqEventObject(fRawMetName,        fRawMet,        true);
  ReqEventObject(fMetMVAName,        fMetMVA,        fMetMVAFromBranch);

  // Initialize the PU histrograms and weights
  // some useful definitions
  if (! fIsData) {
    TString dirFwk("AnaFwkMod");
    TString allEvts("hDAllEvents");
    // get input histo
    TFile *fif = new TFile(fPUInputFileName.Data());
    if (fif->IsOpen() == kFALSE) {
      printf(" WARNING -- missing input pile up file!\n");
    }  
    TDirectory *dirTmp = (TDirectory*) gROOT->FindObject(dirFwk.Data());
    if (dirTmp) {
      fif->cd(dirFwk.Data());
      fPUInput = (TH1D*) dirTmp->Get("hNPUTrue")->Clone();
      if (! fPUInput)
        printf(" WARNING -- no input framework file!\n");      
    }
    // get target histo
    TFile *ftf = new TFile(fPUTargetFileName.Data());
    if (ftf->IsOpen() == kFALSE) {
      printf(" WARNING -- missing target pile up file!\n");
    }  
    fPUTarget = (TH1D*) ftf->Get("pileup")->Clone();
    if (! fPUTarget)
      printf(" WARNING -- no target pile up histogram !\n");      
    // build pile up weight histo
    fPUInput->Rebin(10);
    fPUInput->Scale(1.0/fPUInput->GetSumOfWeights());
    fPUTarget->Scale(1.0/fPUTarget->GetSumOfWeights());
    fPUWeight = new TH1D((*fPUTarget) / (*fPUInput));
  }

  // Create Ntuple Tree
  fOutputFile = TFile::Open(TString::Format("%s_tmp.root",GetName()),"RECREATE");
  fMitDMSTree.CreateTree(0);
  fMitDMSTree.tree_->SetAutoSave(300e9);
  fMitDMSTree.tree_->SetDirectory(fOutputFile);
  AddOutput(fMitDMSTree.tree_);
}

//--------------------------------------------------------------------------------------------------
Float_t DMSTreeWriter::PUWeight(Float_t npu)
{
  if (npu<0)
    return 1.0;
  if (!fPUWeight)
    return 1.0;
  
  return fPUWeight->GetBinContent(fPUWeight->FindFixBin(npu));
}

//--------------------------------------------------------------------------------------------------
void DMSTreeWriter::getGenLevelInfo(MitDMSTree& tree)
{
  // Loop on all stable MC particles
  for (UInt_t i=0; i<fMCParticles->GetEntries(); ++i) {
    const MCParticle *p = fMCParticles->At(i);
    if (p->Status()!=3)
      continue;
    
    // Check if the particle is a Boson
    if (p->Is(MCParticle::kZ) || p->Is(MCParticle::kW)) {
      tree.genV_   = p->Mom();      
      tree.genVid_ = p->AbsPdgId();      

      // Check the daughters of the Boson
      if (p->NDaughters() > 0)
        tree.genVdaughterId_ = p->Daughter(0)->AbsPdgId();
        
      // Special case: top->Wb. Check if the mother is a top
      if (p->Mother()->Is(MCParticle::kTop)) {
        if (p->Mother()->PdgId()==6)      
          tree.topPt_    = p->Pt();
        else 
          tree.topBarPt_ = p->Pt();
      } // end top scope
      
    } // end boson scope 
    // Check if the particle is a high Pt photon
    else if (p->Is(MCParticle::kGamma) && p->Pt() > 100.) {
      tree.genV_   = p->Mom();      
      tree.genVid_ = p->AbsPdgId();            
    }
    else
      continue;
    
  } // end loop on MC Particles

  return;
}
 
//--------------------------------------------------------------------------------------------------
Bool_t DMSTreeWriter::IsHLTMatched(LorentzVector& v,
                                   TriggerObject::ETriggerObject type,
                                   Float_t deltaR,
                                   Float_t minPt,
                                   Float_t maxEta)
{
  float minDr = 999.;
  
  // Loop on the selected trigger object collection 
  // and find a matched object
  for (UInt_t i=0; i<fTrigObj->GetEntries(); ++i) {
    const TriggerObject *to = fTrigObj->At(i);
    if (to->TriggerType() != type)
      continue; 
    // check trigger object pt,eta
    if (to->Pt() < minPt || fabs(to->Eta()) > maxEta)
      continue; 
    // compute the dR
    float thisDr = MathUtils::DeltaR(v, *to);
    if (thisDr < minDr) 
      minDr = thisDr;
    // check if user match condition is met
    if (minDr < deltaR)
      return true;
  } // end loop on trigger objects

  return false;
}


//--------------------------------------------------------------------------------------------------
Bool_t DMSTreeWriter::fjetIsCleaned(LorentzVector& v,
                                    Float_t deltaR)
{
  // Loop on tight muons and photons and discard 
  // jet if overlapping
  for (UInt_t i=0; i<fMuons->GetEntries(); ++i) {
    // Discard loose muons
    if (fMuons->At(i)->ParticleId() == XsIsoParticle::EParticleId::eX)
      continue;
    float thisDr = MathUtils::DeltaR(v, fMuons->At(i)->Mom());
    if (thisDr < deltaR)
      return false;      
  }
  for (UInt_t i=0; i<fPhotons->GetEntries(); ++i) {
    // Discard soft photons
    if (fPhotons->At(i)->Pt() < 150.)
      continue;
    float thisDr = MathUtils::DeltaR(v, fPhotons->At(i)->Mom());
    if (thisDr < deltaR)
      return false;      
  }
 
  return true;      
}

//--------------------------------------------------------------------------------------------------
Int_t DMSTreeWriter::JetPartonMatch(LorentzVector& v,
                                    Float_t deltaR)
{
  float minDr = 999.;
  unsigned int pId = 0;
      
  // Loop on all stable MC particles and perform the matching
  for (UInt_t i=0; i<fMCParticles->GetEntries(); ++i) {
    const MCParticle *p = fMCParticles->At(i);
    if (p->Status()!=3)
      continue;
    // compute the dR    
    float thisDr = MathUtils::DeltaR(v, *p);
    // standard check for any parton    
    if (thisDr < minDr) {
      minDr = thisDr;
      pId = p->AbsPdgId();
    }
    // special check: if quark inside matching cone and mother is W/Z
    // make extra checks and if true return the mother id
    if (thisDr < deltaR && p->AbsPdgId() < 6 
     && (p->Mother()->Is(MCParticle::kZ) || p->Mother()->Is(MCParticle::kW))) {
      // check that both the quarks are inside the matching cone
      float q1Dr = MathUtils::DeltaR(v, *(p->Mother()->Daughter(0)));
      float q2Dr = MathUtils::DeltaR(v, *(p->Mother()->Daughter(1)));
      if (q1Dr < deltaR && q2Dr < deltaR)
        return p->Mother()->AbsPdgId();
    }
  }

  // return matched parton only if minDr less that user dr limit
  if (minDr < deltaR)
    return pId;
  else 
    return 0;  
}

//--------------------------------------------------------------------------------------------------
Float_t DMSTreeWriter::GetFatJetBtag(LorentzVector& v,
                                     Float_t deltaR)
{
  float minDr = 999.;
    
  // Loop on standard jet collection, find the best matched jet 
  // and get the btagging
  for (UInt_t i = 0; i < fJets->GetEntries(); ++i) {
    const Jet *jet = fJets->At(i);
    // compute the dR
    float thisDr = MathUtils::DeltaR(v, *jet);
    if (thisDr < minDr) 
      minDr = thisDr;
    // if the fat jet is matched to standard jet return the btagging
    // of the standard jet
    if (minDr < deltaR)
      return jet->CombinedSecondaryVertexBJetTagsDisc();  
  }
      
  return -1;
}

//--------------------------------------------------------------------------------------------------
void DMSTreeWriter::CorrectMet(const float met, const float metPhi, const LorentzVector& l1, const LorentzVector& l2,
                               float& newMet, float& newMetPhi)
{
  // inputs:  met, metPhi, l1, l2
  // outputs: newMet, newMetPhi
  float newMetX;
  float newMetY;
  newMetX = met*TMath::Cos(metPhi) + l1.Px() + l2.Px();
  newMetY = met*TMath::Sin(metPhi) + l1.Py() + l2.Py();
  
  newMet    = TMath::Sqrt(TMath::Power(newMetX,2) + TMath::Power(newMetY,2));
  newMetPhi = TMath::ATan2(newMetY,newMetX);

  return;
}

//--------------------------------------------------------------------------------------------------
void DMSTreeWriter::CorrectMet(const float met, const float metPhi, const LorentzVector& l1,
                               float& newMet, float& newMetPhi)
{
  // inputs:  met, metPhi, l1
  // outputs: newMet, newMetPhi
  float newMetX;
  float newMetY;
  newMetX = met*TMath::Cos(metPhi) + l1.Px();
  newMetY = met*TMath::Sin(metPhi) + l1.Py();
  
  newMet    = TMath::Sqrt(TMath::Power(newMetX,2) + TMath::Power(newMetY,2));
  newMetPhi = TMath::ATan2(newMetY,newMetX);

  return;
}

//--------------------------------------------------------------------------------------------------
float DMSTreeWriter::GetMt(const LorentzVector& l1, const float met, const float metPhi)
{
  // inputs:  l1, met, metPhi
  // outputs: transverse mass
  double lepPhi = l1.Phi();
  return TMath::Sqrt(2*met*l1.Pt()*(1-TMath::Cos(MathUtils::DeltaPhi(lepPhi,(double)metPhi))));
}

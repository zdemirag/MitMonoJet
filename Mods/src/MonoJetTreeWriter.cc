#include "MitMonoJet/Mods/interface/MonoJetTreeWriter.h"
#include "MitAna/DataTree/interface/PhotonCol.h"
#include "MitAna/DataTree/interface/PFCandidateCol.h"
#include "MitAna/DataTree/interface/GenericParticle.h"
#include "MitAna/DataTree/interface/MCParticleFwd.h"
#include "MitAna/DataTree/interface/PFMet.h"
#include "MitPhysics/Init/interface/ModNames.h"
#include "MitPhysics/Utils/interface/IsolationTools.h"
#include "MitPhysics/Utils/interface/MuonTools.h"
#include "MitPhysics/Utils/interface/MuonIDMVA.h"
#include "MitPhysics/Utils/interface/PhotonTools.h"
#include "MitPhysics/Utils/interface/VertexTools.h"
#include "MitPhysics/Utils/interface/PFMetCorrectionTools.h"
#include "MitAna/DataTree/interface/PFJetCol.h"
#include "MitAna/DataTree/interface/PFTauCol.h"
#include "MitAna/DataTree/interface/GenJetCol.h"
#include "MitPhysics/Mods/interface/MuonIDMod.h"
#include "MitCommon/MathTools/interface/MathUtils.h"
#include "MitAna/DataTree/interface/MuonFwd.h"
#include "MitAna/DataTree/interface/ElectronFwd.h"
#include "MitAna/DataTree/interface/VertexCol.h"
#include "MitAna/TreeMod/interface/BaseMod.h" 
#include "TDataMember.h"
#include "TFile.h"
#include <TNtuple.h>
#include <TRandom3.h>
#include <TSystem.h>

#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"

using namespace mithep;

ClassImp(mithep::MonoJetTreeWriter)

//--------------------------------------------------------------------------------------------------
MonoJetTreeWriter::MonoJetTreeWriter(const char *name, const char *title) : 
  // Base Module...
  BaseMod                 (name,title),
  
  fMetName                ("PFMet"),
  fPhotonsName            (Names::gkPhotonBrn),
  fElectronsName          (Names::gkElectronBrn),
  fMuonsName              (Names::gkMuonBrn),
  fTausName               (Names::gkPFTauBrn),
  fJetsName               (Names::gkPFJetBrn),
  fLeptonsName            (ModNames::gkMergedLeptonsName),
  fPFCandidatesName       (Names::gkPFCandidatesBrn),
  fVertexName             (ModNames::gkGoodVertexesName),

  fSuperClustersName      ("PFSuperClusters"),
  fTracksName             (Names::gkTrackBrn),
  fPVName                 (Names::gkPVBeamSpotBrn),
  fPileUpDenName          (Names::gkPileupEnergyDensityBrn),
  fPileUpName             (Names::gkPileupInfoBrn),  
  fBeamspotName           (Names::gkBeamSpotBrn),
  fMCEvInfoName           (Names::gkMCEvtInfoBrn),

  fMCPartName(Names::gkMCPartBrn),

  fIsData                 (false),
  fTriggerObjectsName     ("MyHltPhotObjs"),
  fMetFromBranch          (kTRUE),  
  fPhotonsFromBranch      (kTRUE),  
  fElectronsFromBranch    (kTRUE),  
  fMuonsFromBranch        (kTRUE),  
  fTausFromBranch         (kTRUE),  
  fJetsFromBranch         (kTRUE),
  fPFCandidatesFromBranch  (kTRUE),
  fPVFromBranch           (kTRUE),
  fQGTaggerCHS            (kTRUE),


  // ----------------------------------------
  fRawMet                 (0),
  fMet                    (0),
  fPhotons                (0),
  fElectrons              (0),
  fMuons                  (0),
  fPFTaus                 (0),
  fJets                   (0),
  fJetCorrector           (0),
  fJetUncertainties       (0),
  fTrigObj                (0),
  fPFCandidates           (0),
  fTracks                 (0),
  fPV                     (0),
  fBeamspot               (0),
  fMCEventInfo            (0),
  fPileUp                 (0),
  fPileUpDen              (0),
  fSuperClusters          (0),
  fParticles              (0),
  fEvtSelData             (0),

  fDecay(0),
  fOutputFile(0),
  fFillNtupleType(0),
  fNEventsSelected(0)

{
  // WARNING, defining the object here invalidates the call of the setter for the CHS flag
  qgTagger = new QGTagger(fQGTaggerCHS);

  // Constructor
}

MonoJetTreeWriter::~MonoJetTreeWriter()
{
  // Destructor
  fOutputFile->Close();
}

//--------------------------------------------------------------------------------------------------
void MonoJetTreeWriter::SlaveTerminate()
{
  fOutputFile->WriteTObject(fMitGPTree.tree_,fMitGPTree.tree_->GetName());
  cout << "Processed events on MonoJetTreeWriter: " << fNEventsSelected << endl;
  delete fJetCorrector;
  delete fJetUncertainties;
  delete fMVAMet;
}


//--------------------------------------------------------------------------------------------------
void MonoJetTreeWriter::Process()
{

  // ------------------------------------------------------------  
  // Process entries of the tree. 
  LoadEventObject("PFMet",            fRawMet,        true);
  LoadEventObject(fMetName,           fMet,           fMetFromBranch);
  LoadEventObject(fPhotonsName,       fPhotons,       fPhotonsFromBranch); 
  LoadEventObject(fElectronsName,     fElectrons,     fElectronsFromBranch);
  LoadEventObject(fMuonsName,         fMuons,         fMuonsFromBranch);
  LoadEventObject(fTausName,          fPFTaus,        fTausFromBranch);
  LoadEventObject(fJetsName,          fJets,          fJetsFromBranch);
  LoadEventObject(fPFCandidatesName,   fPFCandidates, fPFCandidatesFromBranch);

  LoadEventObject(fPVName,            fPV,            fPVFromBranch);    
  LoadEventObject(fBeamspotName,      fBeamspot);
  
  LoadEventObject(fSuperClustersName, fSuperClusters);
  LoadEventObject(fTracksName,        fTracks,        true);

  LoadEventObject(fPileUpDenName,     fPileUpDen,     true);
  if (!fIsData) {
    ReqBranch(fPileUpName,            fPileUp);
    LoadEventObject(fMCPartName, fParticles);
  }

  ParticleOArr *leptons = GetObjThisEvt<ParticleOArr>(ModNames::gkMergedLeptonsName);
  //const MuonCol *muons = GetObjThisEvt<MuonCol>("HggLeptonTagMuons"); //This should be identical to MuonIDMod->GetOutputName() in run macro
  const MuonCol *muons = GetObjThisEvt<MuonCol>(fMuonsName);
  const VertexCol *vertices = GetObjThisEvt<VertexOArr>(fVertexName);
  LoadEventObject("EvtSelData",       fEvtSelData,    true);

  fNEventsSelected++;

  fMitGPTree.InitVariables();
      
  // Pileup related stuff
  if( !fIsData ) {
    LoadBranch(fPileUpName);
    for (UInt_t i=0; i<fPileUp->GetEntries(); ++i) {
      const PileupInfo *puinfo = fPileUp->At(i);
      if (puinfo->GetBunchCrossing() ==  0) fMitGPTree.npu_	    = puinfo->GetPU_NumMean();
      if (puinfo->GetBunchCrossing() ==  1) fMitGPTree.npuPlusOne_  = puinfo->GetPU_NumInteractions();
      if (puinfo->GetBunchCrossing() == -1) fMitGPTree.npuMinusOne_ = puinfo->GetPU_NumInteractions();
    }
  }

  // Trigger stuff

  //cout<<"Event Number: "<<fNEventsSelected<<endl;
  fMitGPTree.trigger_  = 0; 
  fTrigObj = GetHLTObjects(fTriggerObjectsName);
  int n_trigjets=0;
  if (!fTrigObj) printf("MonoJetTreeWriter::TriggerObjectCol not found\n");
  else {
    for(UInt_t i=0;i<fTrigObj->GetEntries();++i) {
      const TriggerObject *to = fTrigObj->At(i);
      TString trName = to->TrigName();
     
      /*if(trName.Contains("HLT_DiPFJet40_PFMETnoMu65_MJJ800VBF_AllJets_v"))
	//if(trName.Contains("HLT_DiPFJet40_PFMETnoMu65_MJJ600VBF_LeadingJets_v"))
	{
	  if(to->Type()==85)
	    {
	      cout<<"Trigger Object Id: "<<to->Type()<<endl;
	      cout<<"Trigger Object pT: "<<to->Pt()<<endl;
	      cout<<"Trigger Object eta: "<<to->Eta()<<endl;
	      cout<<"Trigger Object phi: "<<to->Phi()<<endl;
	      cout<<"Trigger Object Name: "<<to->TrigName()<<endl<<endl<<endl;
	      n_trigjets++;

	      //cout<<"Number of trigger Jets: "<<n_trigjets<<endl;
	    }
	    }*/

      // default MonoJet
      if(trName.Contains("MonoCentralPFJet80_PFMETnoMu"))
	fMitGPTree.trigger_ |= 1 << 0;
      if(trName.Contains("HLT_MET120_HBHENoiseCleaned_v") )
      fMitGPTree.trigger_ |= 1 << 1;
      // default VBF
      if(trName.Contains("HLT_DiPFJet40_PFMETnoMu65_MJJ800VBF_AllJets_v")) 
	fMitGPTree.trigger_ |= 1 << 2;
      // parked VBF, B-D. FIXME
      if(trName.Contains("HLT_DiJet30_MJJ700_AllJets_DEta3p5_VBF_v1") or 
	 trName.Contains("HLT_DiJet35_MJJ700_AllJets_DEta3p5_VBF_v5"))	
	fMitGPTree.trigger_ |= 1 << 3;
      // default single muon
      if(trName.Contains("HLT_IsoMu24_eta2p1_v")) 
      fMitGPTree.trigger_ |= 1 << 4;
    }
  }

  //cout<<"Number of trigger Jets: "<<n_trigjets<<endl;

  fMitGPTree.run_   = GetEventHeader()->RunNum();
  fMitGPTree.lumi_  = GetEventHeader()->LumiSec();
  fMitGPTree.event_ = GetEventHeader()->EvtNum();
  fMitGPTree.nvtx_  = fPV->GetEntries();
  fMitGPTree.scale1fb_ = 1000.0;
  fMitGPTree.cuts_ = MitGPTree::DiLepton;
  
  if(fDecay == 0) fMitGPTree.dstype_ = MitGPTree::data;
  else            fMitGPTree.dstype_ = MitGPTree::other;

  fMitGPTree.metRaw_    = fRawMet->At(0)->Pt();
  fMitGPTree.metRawPhi_ = fRawMet->At(0)->Phi();
  fMitGPTree.metRawCorZ_    = fRawMet->At(0)->Pt();
  fMitGPTree.metRawCorZPhi_ = fRawMet->At(0)->Phi();
  fMitGPTree.metRawCorW_    = fRawMet->At(0)->Pt();
  fMitGPTree.metRawCorWPhi_ = fRawMet->At(0)->Phi();
  fMitGPTree.met_       = fMet->At(0)->Pt();
  fMitGPTree.metPhi_    = fMet->At(0)->Phi();
  fMitGPTree.metCorZ_   = fMet->At(0)->Pt();  //MetCor default value as Met
  fMitGPTree.metCorZPhi_= fMet->At(0)->Phi(); //MetCor default value as Met
  fMitGPTree.metCorW_   = fMet->At(0)->Pt();  //MetCor default value as Met
  fMitGPTree.metCorWPhi_= fMet->At(0)->Phi(); //MetCor default value as Met
  fMitGPTree.sumEt_     = fMet->At(0)->SumEt();
  fMitGPTree.metSig_    = fRawMet->At(0)->PFMetSig(); //Use the RawPFMet for the sig. calculation


  // TAUs
  fMitGPTree.ntaus_ = fPFTaus->GetEntries();
  if (fPFTaus->GetEntries() >= 1) {
    const PFTau *tau = fPFTaus->At(0);
    fMitGPTree.tau1_ = tau->Mom();
    if (fPFTaus->GetEntries() >= 2) {
      tau = fPFTaus->At(1);
      fMitGPTree.tau2_ = tau->Mom();
    }
  }

  // LEPTONS
  unsigned int muoncount = 0;
  const Muon *mu;
  fMitGPTree.nlep_ = leptons->GetEntries();
  if (leptons->GetEntries() >= 1) {
    const Particle *lep = leptons->At(0);
    if(muons->GetEntries() > muoncount){ mu = muons->At(muoncount);}
    
    fMitGPTree.lep1_ = lep->Mom();
    if     (lep->ObjType() == kMuon    ){
    muoncount++;
    fMitGPTree.lid1_ = 13;
    if(((mu->HasGlobalTrk() && mu->GlobalTrk()->Chi2()/mu->GlobalTrk()->Ndof() < 10 && (mu->NSegments() > 1 || mu->NMatches() > 1) && mu->NValidHits() > 0) 
	|| mu->IsTrackerMuon()) && (mu->BestTrk() != 0 && mu->BestTrk()->NHits() > 10 &&
				    (mu->NSegments() > 1 || mu->NMatches() > 1) && mu->BestTrk()->NPixelHits() > 0))
      { fMitGPTree.lep1IsTightMuon_ = 1;}
    }
    else if(lep->ObjType() == kElectron) fMitGPTree.lid1_ = 11;
    else                                 assert(0);
    if(lep->Charge() < 0) fMitGPTree.lid1_ = -1 * fMitGPTree.lid1_;
    // If the event contains more than 1 leptons correct the MET using the highest pt one
    float met_new_x = fMitGPTree.met_*TMath::Cos(fMitGPTree.metPhi_) + fMitGPTree.lep1_.Px();
    float met_new_y = fMitGPTree.met_*TMath::Sin(fMitGPTree.metPhi_) + fMitGPTree.lep1_.Py();

    fMitGPTree.metCorW_    = TMath::Sqrt(TMath::Power(met_new_x,2) + TMath::Power(met_new_y,2));
    fMitGPTree.metCorWPhi_ = TMath::ATan2(met_new_y,met_new_x);

    float metRaw_new_x = fMitGPTree.metRaw_*TMath::Cos(fMitGPTree.metRawPhi_) + fMitGPTree.lep1_.Px();
    float metRaw_new_y = fMitGPTree.metRaw_*TMath::Sin(fMitGPTree.metRawPhi_) + fMitGPTree.lep1_.Py();

    fMitGPTree.metRawCorW_    = TMath::Sqrt(TMath::Power(metRaw_new_x,2) + TMath::Power(metRaw_new_y,2));
    fMitGPTree.metRawCorWPhi_ = TMath::ATan2(metRaw_new_y,metRaw_new_x);
  }
  if (leptons->GetEntries() >= 2) {
    const Particle *lep = leptons->At(1);
    if(muons->GetEntries() > muoncount){ mu = muons->At(muoncount);}
    
    fMitGPTree.lep2_ = lep->Mom();
    if     (lep->ObjType() == kMuon    ){
    muoncount++;
    fMitGPTree.lid2_ = 13;
    if(((mu->HasGlobalTrk() && mu->GlobalTrk()->Chi2()/mu->GlobalTrk()->Ndof() < 10 &&
	       (mu->NSegments() > 1 || mu->NMatches() > 1) && mu->NValidHits() > 0) ||
	       mu->IsTrackerMuon()) && (mu->BestTrk() != 0 && mu->BestTrk()->NHits() > 10 &&
		(mu->NSegments() > 1 || mu->NMatches() > 1) && mu->BestTrk()->NPixelHits() > 0)){ fMitGPTree.lep2IsTightMuon_ = 1;}
    }
    else if(lep->ObjType() == kElectron) fMitGPTree.lid2_ = 11;
    else                                 assert(0);
    if(lep->Charge() < 0) fMitGPTree.lid2_ = -1 * fMitGPTree.lid2_;
    // If the event contains more than 2 leptons correct the MET using the 2 highest pt ones
    float met_new_x = fMitGPTree.met_*TMath::Cos(fMitGPTree.metPhi_) + fMitGPTree.lep1_.Px() + fMitGPTree.lep2_.Px();
    float met_new_y = fMitGPTree.met_*TMath::Sin(fMitGPTree.metPhi_) + fMitGPTree.lep1_.Py() + fMitGPTree.lep2_.Py();

    fMitGPTree.metCorZ_    = TMath::Sqrt(TMath::Power(met_new_x,2) + TMath::Power(met_new_y,2));
    fMitGPTree.metCorZPhi_ = TMath::ATan2(met_new_y,met_new_x);

    float metRaw_new_x = fMitGPTree.metRaw_*TMath::Cos(fMitGPTree.metRawPhi_) + fMitGPTree.lep1_.Px() + fMitGPTree.lep2_.Px();
    float metRaw_new_y = fMitGPTree.metRaw_*TMath::Sin(fMitGPTree.metRawPhi_) + fMitGPTree.lep1_.Py() + fMitGPTree.lep2_.Py();

    fMitGPTree.metRawCorZ_    = TMath::Sqrt(TMath::Power(metRaw_new_x,2) + TMath::Power(metRaw_new_y,2));
    fMitGPTree.metRawCorZPhi_ = TMath::ATan2(metRaw_new_y,metRaw_new_x);
  }
  if (leptons->GetEntries() >= 3) {
    const Particle *lep = leptons->At(2);
    if(muons->GetEntries() > muoncount){ mu = muons->At(muoncount);}
    
    fMitGPTree.lep3_ = lep->Mom();
    if     (lep->ObjType() == kMuon    ){
    muoncount++;
    fMitGPTree.lid3_ = 13;
    if(((mu->HasGlobalTrk() && mu->GlobalTrk()->Chi2()/mu->GlobalTrk()->Ndof() < 10 &&
	       (mu->NSegments() > 1 || mu->NMatches() > 1) && mu->NValidHits() > 0) ||
	       mu->IsTrackerMuon()) && (mu->BestTrk() != 0 && mu->BestTrk()->NHits() > 10 &&
		(mu->NSegments() > 1 || mu->NMatches() > 1) && mu->BestTrk()->NPixelHits() > 0)){ fMitGPTree.lep3IsTightMuon_ = 1;}
    }
    else if(lep->ObjType() == kElectron) fMitGPTree.lid3_ = 11;
    else                                 assert(0);
    if(lep->Charge() < 0) fMitGPTree.lid3_ = -1 * fMitGPTree.lid3_;
  }

  //PHOTONS  
  fMitGPTree.nphotons_ = fPhotons->GetEntries();
  if(fPhotons->GetEntries() >= 1) {
    const Photon *photon = fPhotons->At(0);
    fMitGPTree.pho1_			  = photon->Mom();
  }
//   if(fPhotons->GetEntries() >= 2) {
//     const Photon *photon = fPhotons->At(1);
//     fMitGPTree.pho2_			  = photon->Mom();
//   }
//   if(fPhotons->GetEntries() >= 3) {
//     const Photon *photon = fPhotons->At(2);
//     fMitGPTree.pho3_			  = photon->Mom();
//   }
//   if(fPhotons->GetEntries() >= 4) {
//     const Photon *photon = fPhotons->At(3);
//     fMitGPTree.pho4_			  = photon->Mom();
//   }

  //JETS
  bool jet1HLTmatch=false;
  double trigjet1_Px=0.;
  double trigjet1_Py=0.;
  double trigjet1_Pz=0.;
  double trigjet1_E=0.;
  double trigjet1_Eta=0.;

  bool jet2HLTmatch=false;
  double trigjet2_Px=0.;
  double trigjet2_Py=0.;
  double trigjet2_Pz=0.;
  double trigjet2_E=0.;
  double trigjet2_Eta=0.;

  PFJetOArr *PFJets = new PFJetOArr;
  PFJets->SetOwner(kTRUE);
  //PFJets->SetName(fPFJetsName);

  for (UInt_t i=0; i<fJets->GetEntries(); ++i) {
    const Jet *inJet = fJets->At(i);
  
    //copy input jet, using special function to copy full derived class
    Jet *jet = inJet->MakeCopy();
    PFJets->AddOwned(dynamic_cast<PFJet*>(jet));

    //cache uncorrected momentum
    const FourVectorM rawMom = jet->RawMom();

    //compute correction factors
    fJetCorrector->setJetEta(rawMom.Eta());
    fJetCorrector->setJetPt(rawMom.Pt());
    fJetCorrector->setJetPhi(rawMom.Phi());
    fJetCorrector->setJetE(rawMom.E());
    
    fJetCorrector->setRho(fPileUpDen->At(0)->RhoRandom());
    fJetCorrector->setJetA(jet->JetArea());

    fJetCorrector->setJetEMF(-99.0);

  }  

  fMitGPTree.njets_ = fJets->GetEntries();
  fMitGPTree.noiseCleaning_ = 0;
  //qgTagger->SetRhoIso(fPileUpDen->At(0)->RhoKt6PFJetsCentralChargedPileUp()); // is it this one?
  qgTagger->SetRhoIso(fPileUpDen->At(0)->RhoRandomLowEta()); // is it this one?
  if (fJets->GetEntries() >= 1) {
    const PFJet *jet = dynamic_cast<const PFJet*>(fJets->At(0));
    fMitGPTree.jet1_     = jet->Mom();
    fJetUncertainties->setJetPt(jet->Pt());
    fJetUncertainties->setJetEta(jet->Eta());
    fMitGPTree.jet1Unc_ = fJetUncertainties ->getUncertainty(true);
    fMitGPTree.jet1CHF_  = jet->ChargedHadronEnergy()/jet->RawMom().E();
    fMitGPTree.jet1NHF_  = jet->NeutralHadronEnergy()/jet->RawMom().E();
    fMitGPTree.jet1NEMF_  = jet->NeutralEmEnergy()/jet->RawMom().E();
    fMitGPTree.noiseCleaning_ |= int(fMitGPTree.jet1CHF_>0.2) << 0;
    fMitGPTree.noiseCleaning_ |= int(fMitGPTree.jet1NHF_>0.7) << 1;
    fMitGPTree.noiseCleaning_ |= int(fMitGPTree.jet1NEMF_>0.7) << 2;
    fMitGPTree.jet1Btag_ = jet->CombinedSecondaryVertexBJetTagsDisc();
    qgTagger->CalculateVariables(jet, vertices);
    fMitGPTree.jet1QGtag_ = qgTagger->QGValue();

    // variables used for the QG retraining
    fMitGPTree.jet1QGRho_   = fPileUpDen->At(0)->RhoRandomLowEta();
    fMitGPTree.jet1QGPtD_   = qgTagger->GetPtD();
    fMitGPTree.jet1QGAxis1_ = qgTagger->GetAxis1();
    fMitGPTree.jet1QGAxis2_ = qgTagger->GetAxis2();
    fMitGPTree.jet1QGMult_  = qgTagger->GetMult();
    
    // matching
    if (!fIsData) {
      double minPartonicDR = 0.8;
      unsigned int partonId = 0;
      for (UInt_t i=0; i<fParticles->GetEntries(); ++i) {
        const MCParticle *p = fParticles->At(i);
        if (p->Status()!=3) continue;
        if (p->AbsPdgId()>5 and p->AbsPdgId()!=21) continue; //1, 2, 3, 4, 5, 21 
        if (MathUtils::DeltaR(*p,*jet)< minPartonicDR){
          minPartonicDR = MathUtils::DeltaR(*p,*jet);
          partonId = p->AbsPdgId();
        }
      }
      fMitGPTree.jet1PartonId_ = partonId;
    }

    // trigger matching
    for(unsigned int i=0; i<fTrigObj->GetEntries(); ++i){
      const TriggerObject *trigobj=fTrigObj->At(i);
      TString trigName = trigobj->TrigName();
      if(trigobj->IsHLT() ){
	if(trigName.Contains("MonoCentralPFJet80_PFMETnoMu"))
	  {
	    bool match = true;
	    if(trigobj->Pt()<80) match=false;
	    if(trigobj->Type()!=85) match=false;
	    if(MathUtils::DeltaR(trigobj->Phi(), trigobj->Eta(), jet->Phi(), jet->Eta())>0.5) match=false;
	    if(match)
	      {
		fMitGPTree.HLTmatch_ |= 1 << 0;
	      }
	  }
	if(trigName.Contains("HLT_DiPFJet40_PFMETnoMu65_MJJ800VBF_AllJets_v"))
	  {
	    bool match = true;
	   
	    if(trigobj->Pt()<40) match=false;
	    if(trigobj->Type()!=85) match=false;
	    if(MathUtils::DeltaR(trigobj->Phi(), trigobj->Eta(), jet->Phi(), jet->Eta())>0.5) match=false;
	    if(match)
	      {
		jet1HLTmatch=true;
		trigjet1_Px=trigobj->Px();
		trigjet1_Py=trigobj->Py();
		trigjet1_Pz=trigobj->Pz();
		trigjet1_E=trigobj->E();
		trigjet1_Eta=trigobj->Eta();
	      }
	  }
      } 
    }
  }

  if (fJets->GetEntries() >= 2) {
    const PFJet *jet = dynamic_cast<const PFJet*>(fJets->At(1));
    fMitGPTree.jet2_     = jet->Mom();
    fJetUncertainties->setJetPt(jet->Pt());
    fJetUncertainties->setJetEta(jet->Eta());
    fMitGPTree.jet2Unc_ = fJetUncertainties ->getUncertainty(true);
    fMitGPTree.jet2CHF_  = jet->ChargedHadronEnergy()/jet->RawMom().E();
    fMitGPTree.jet2NHF_  = jet->NeutralHadronEnergy()/jet->RawMom().E();
    fMitGPTree.jet2NEMF_  = jet->NeutralEmEnergy()/jet->RawMom().E();
    fMitGPTree.noiseCleaning_ |= int(fMitGPTree.jet2CHF_>0.2) << 3;
    fMitGPTree.noiseCleaning_ |= int(fMitGPTree.jet2NHF_>0.7) << 4;
    fMitGPTree.noiseCleaning_ |= int(fMitGPTree.jet2NEMF_>0.7) << 5;
    fMitGPTree.jet2Btag_ = jet->CombinedSecondaryVertexBJetTagsDisc();
    qgTagger->CalculateVariables(jet, vertices);
    fMitGPTree.jet2QGtag_ = qgTagger->QGValue();


    // trigger matching
    for(unsigned int i=0; i<fTrigObj->GetEntries(); ++i){
      const TriggerObject *trigobj=fTrigObj->At(i);
      TString trigName = trigobj->TrigName();
      if(trigobj->IsHLT() ){
	if(trigName.Contains("HLT_DiPFJet40_PFMETnoMu65_MJJ800VBF_AllJets_v"))
	  {
	    bool match = true;
	    if(trigobj->Pt()<40) match=false;
	    if(trigobj->Type()!=85) match=false;
	    if(MathUtils::DeltaR(trigobj->Phi(), trigobj->Eta(), jet->Phi(), jet->Eta())>0.5) match=false;
	    if(match)
	      {
		jet2HLTmatch=true;
		trigjet2_Px=trigobj->Px();
		trigjet2_Py=trigobj->Py();
		trigjet2_Pz=trigobj->Pz();
		trigjet2_E=trigobj->E();
		trigjet2_Eta=trigobj->Eta();
	      }
	  }
      } 
    } 
  }

  if(jet1HLTmatch&&jet2HLTmatch)
    {
      double dijetmass=sqrt(pow(trigjet1_E+trigjet2_E,2)-(pow(trigjet1_Px+trigjet2_Px,2)+pow(trigjet1_Py+trigjet2_Py,2)+pow(trigjet1_Pz+trigjet2_Pz,2)));
      double deltaeta=fabs(trigjet1_Eta-trigjet2_Eta);
      if((dijetmass>800)&&(deltaeta>3.5)&&(trigjet1_Eta*trigjet2_Eta<0))
	{
	  fMitGPTree.HLTmatch_ |= 1 << 2;
	}
    }

  if (fJets->GetEntries() >= 3) {
    const PFJet *jet = dynamic_cast<const PFJet*>(fJets->At(2));
    fMitGPTree.jet3_     = jet->Mom();
    fMitGPTree.jet3CHF_  = jet->ChargedHadronEnergy()/jet->RawMom().E();
    fMitGPTree.jet3NHF_  = jet->NeutralHadronEnergy()/jet->RawMom().E();
    fMitGPTree.jet3NEMF_  = jet->NeutralEmEnergy()/jet->RawMom().E();
    fMitGPTree.jet3Btag_ = jet->CombinedSecondaryVertexBJetTagsDisc();
  }
  if (fJets->GetEntries() >= 4) {
    const Jet *jet = fJets->At(3);
    fMitGPTree.jet4_     = jet->Mom();
    fMitGPTree.jet4Btag_ = jet->CombinedSecondaryVertexBJetTagsDisc();
  }

          
  //TRACKS
  fMitGPTree.ntracks_ = 0;
  for(unsigned int i = 0; i < fTracks->GetEntries(); i++) {
    const mithep::Track* pTrack = fTracks->At(i);
    if(pTrack->Pt() <= 15) continue;
    Bool_t isLepton = kFALSE;
    for (unsigned int arrayIndex=0; arrayIndex<leptons->GetEntries(); arrayIndex++ ) {
       const Particle *lep = leptons->At(arrayIndex);
      if(MathUtils::DeltaR(pTrack->Mom(), lep->Mom()) < 0.05) {
        isLepton = kTRUE;
        break;
      }
    }
    if(isLepton == kTRUE) continue;
    GenericParticle *p = new GenericParticle(pTrack->Px(), pTrack->Py(), pTrack->Pz(), 
                                             pTrack->P(), pTrack->Charge());
    if(fMitGPTree.ntracks_ == 0) fMitGPTree.track1_ = p->Mom();
    if(fMitGPTree.ntracks_ == 1) fMitGPTree.track2_ = p->Mom();
    if(fMitGPTree.ntracks_ == 2) fMitGPTree.track3_ = p->Mom();
    delete p;
    fMitGPTree.ntracks_++;
  }

  Met mvaMet = fMVAMet->GetMet(fMuons,fElectrons,fPFTaus,fPFCandidates,PFJets,0,fPV,fRawMet,fJetCorrector,fPileUpDen);
  TMatrixD* MVACov = fMVAMet->GetMetCovariance();

  fMitGPTree.mvamet_ = mvaMet.Pt();
  fMitGPTree.mvametPhi_ = mvaMet.Phi();
  fMitGPTree.mvametCorZ_ = mvaMet.Pt();
  fMitGPTree.mvametCorZPhi_ = mvaMet.Phi();
  fMitGPTree.mvametCorW_ = mvaMet.Pt();
  fMitGPTree.mvametCorWPhi_ = mvaMet.Phi();
  fMitGPTree.mvaCov00_ = (*MVACov)(0,0);
  fMitGPTree.mvaCov10_ = (*MVACov)(1,0);
  fMitGPTree.mvaCov01_ = (*MVACov)(0,1);
  fMitGPTree.mvaCov11_ = (*MVACov)(1,1);

  if (leptons->GetEntries() >= 1) {
    float mvamet_new_x = fMitGPTree.mvamet_*TMath::Cos(fMitGPTree.mvametPhi_) + fMitGPTree.lep1_.Px();
    float mvamet_new_y = fMitGPTree.mvamet_*TMath::Sin(fMitGPTree.mvametPhi_) + fMitGPTree.lep1_.Py();

    fMitGPTree.mvametCorW_    = TMath::Sqrt(TMath::Power(mvamet_new_x,2) + TMath::Power(mvamet_new_y,2));
    fMitGPTree.mvametCorWPhi_ = TMath::ATan2(mvamet_new_y,mvamet_new_x);
  }
  if (leptons->GetEntries() >= 2) {
    float mvamet_new_x = fMitGPTree.mvamet_*TMath::Cos(fMitGPTree.mvametPhi_) + fMitGPTree.lep1_.Px() + fMitGPTree.lep2_.Px();
    float mvamet_new_y = fMitGPTree.mvamet_*TMath::Sin(fMitGPTree.mvametPhi_) + fMitGPTree.lep1_.Py() + fMitGPTree.lep2_.Py();

    fMitGPTree.mvametCorZ_    = TMath::Sqrt(TMath::Power(mvamet_new_x,2) + TMath::Power(mvamet_new_y,2));
    fMitGPTree.mvametCorZPhi_ = TMath::ATan2(mvamet_new_y,mvamet_new_x);
  }

  Double_t Q	     = 0.0;
  Int_t    id1       = 0;
  Double_t x1	     = 0.0;
  Double_t pdf1      = 0.0;
  Int_t    id2       = 0;
  Double_t x2	     = 0.0;
  Double_t pdf2      = 0.0;
  Int_t    processId = 0;
  if(fIsData == kFALSE){
     LoadBranch(fMCEvInfoName);
     Q         = fMCEventInfo->Scale();
     id1       = fMCEventInfo->Id1();
     x1        = fMCEventInfo->X1();
     pdf1      = fMCEventInfo->Pdf1();
     id2       = fMCEventInfo->Id2();
     x2        = fMCEventInfo->X2();
     pdf2      = fMCEventInfo->Pdf2();
     processId = fMCEventInfo->ProcessId();
  }
  fMitGPTree.Q_         = Q;
  fMitGPTree.id1_       = id1;
  fMitGPTree.x1_        = x1;
  fMitGPTree.pdf1_      = pdf1;
  fMitGPTree.id2_       = id2;
  fMitGPTree.x2_        = x2;
  fMitGPTree.pdf2_      = pdf2;
  fMitGPTree.processId_ = processId;
  fMitGPTree.metFiltersWord_ = fEvtSelData->metFiltersWord();
  fMitGPTree.tree_->Fill();
  
  return;
}

//--------------------------------------------------------------------------------------------------
void MonoJetTreeWriter::SlaveBegin()
{
  // Run startup code on the computer (slave) doing the actual analysis. Here,
  // we just request the photon collection branch.  
  ReqEventObject("PFMet",            fRawMet,         true);
  ReqEventObject(fMetName,           fMet,            fMetFromBranch);
  ReqEventObject(fPhotonsName,       fPhotons,        fPhotonsFromBranch); 
  ReqEventObject(fElectronsName,     fElectrons,      fElectronsFromBranch);
  ReqEventObject(fMuonsName,         fMuons,          fMuonsFromBranch);
  ReqEventObject(fTausName,          fPFTaus,         fTausFromBranch);
  ReqEventObject(fJetsName,          fJets,           fJetsFromBranch);
  ReqEventObject(fPFCandidatesName,  fPFCandidates,   fPFCandidatesFromBranch);

  ReqEventObject(fSuperClustersName,  fSuperClusters, true);
  ReqEventObject(fTracksName,         fTracks,        true);

  ReqEventObject(fPVName,             fPV,            fPVFromBranch);    
  ReqEventObject(fBeamspotName,       fBeamspot,      true);

  ReqEventObject(fPileUpDenName,   fPileUpDen,    true);
  if (!fIsData) { 
    ReqBranch(fPileUpName,         fPileUp);
    ReqBranch(fMCEvInfoName,       fMCEventInfo);
    ReqEventObject(fMCPartName, fParticles, kTRUE);
  }
  ReqEventObject("EvtSelData",        fEvtSelData,    true);

  if (fIsData )
    {
      fCorrectionFiles.push_back(std::string((gSystem->Getenv("CMSSW_BASE") + TString("/src/MitPhysics/data/Summer13_V1_DATA_L1FastJet_AK5PF.txt")).Data()));
      fCorrectionFiles.push_back(std::string((gSystem->Getenv("CMSSW_BASE") + TString("/src/MitPhysics/data/Summer13_V1_DATA_L2Relative_AK5PF.txt")).Data()));
      fCorrectionFiles.push_back(std::string((gSystem->Getenv("CMSSW_BASE") + TString("/src/MitPhysics/data/Summer13_V1_DATA_L3Absolute_AK5PF.txt")).Data()));
      fCorrectionFiles.push_back(std::string((gSystem->Getenv("CMSSW_BASE") + TString("/src/MitPhysics/data/Summer13_V1_DATA_L2L3Residual_AK5PF.txt")).Data()));
    }
  else
    {
      fCorrectionFiles.push_back(std::string((gSystem->Getenv("CMSSW_BASE") + TString("/src/MitPhysics/data/Summer13_V1_MC_L1FastJet_AK5PF.txt")).Data()));
      fCorrectionFiles.push_back(std::string((gSystem->Getenv("CMSSW_BASE") + TString("/src/MitPhysics/data/Summer13_V1_MC_L2Relative_AK5PF.txt")).Data()));
      fCorrectionFiles.push_back(std::string((gSystem->Getenv("CMSSW_BASE") + TString("/src/MitPhysics/data/Summer13_V1_MC_L3Absolute_AK5PF.txt")).Data()));
    }

  //fill JetCorrectorParameters from files
  std::vector<JetCorrectorParameters> correctionParameters;
  for (std::vector<std::string>::const_iterator it = fCorrectionFiles.begin(); it!=fCorrectionFiles.end(); ++it) {
    correctionParameters.push_back(JetCorrectorParameters(*it));
  }
  
  //initialize jet corrector class
  fJetCorrector = new FactorizedJetCorrector(correctionParameters);

  std::string jetCorrectorParams;
  if (fIsData ) jetCorrectorParams = std::string(TString::Format("%s/src/MitPhysics/data/Summer13_V1_DATA_Uncertainty_AK5PF.txt", getenv("CMSSW_BASE")));
  else jetCorrectorParams = std::string(TString::Format("%s/src/MitPhysics/data/Summer13_V1_MC_Uncertainty_AK5PF.txt", getenv("CMSSW_BASE")));
  JetCorrectorParameters param(jetCorrectorParams);
  fJetUncertainties = new JetCorrectionUncertainty(param);

  fMVAMet     = new MVAMet();

  fMVAMet->Initialize(TString(getenv("CMSSW_BASE")+string("/src/MitPhysics/data/TMVAClassificationCategory_JetID_MET_53X_Dec2012.weights.xml")),
		      TString(getenv("CMSSW_BASE")+string("/src/MitPhysics/data/TMVAClassificationCategory_JetID_MET_53X_Dec2012.weights.xml")),
		      TString(getenv("CMSSW_BASE")+string("/src/MitPhysics/Utils/python/JetIdParams_cfi.py")),
		      TString(getenv("CMSSW_BASE")+string("/src/MitPhysics/data/gbrmet_53_Dec2012.root")),
		      TString(getenv("CMSSW_BASE")+string("/src/MitPhysics/data/gbrmetphi_53_Dec2012.root")),
		      TString(getenv("CMSSW_BASE")+string("/src/MitPhysics/data/gbru1cov_53_Dec2012.root")),
		      TString(getenv("CMSSW_BASE")+string("/src/MitPhysics/data/gbru2cov_53_Dec2012.root")),JetIDMVA::k53MET);

  //***********************************************************************************************
  //Create Smurf Ntuple Tree
  //***********************************************************************************************
  fOutputFile = TFile::Open(TString::Format("%s_tmp.root",GetName()),"RECREATE");
  fMitGPTree.CreateTree(fFillNtupleType);
  fMitGPTree.tree_->SetAutoSave(300e9);
  fMitGPTree.tree_->SetDirectory(fOutputFile);
  AddOutput(fMitGPTree.tree_);

  
}

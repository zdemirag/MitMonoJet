// $Id $

#include <iostream>
#include <sstream>
#include <TH1D.h>
#include <TH2D.h>
#include <TMath.h>

#include "MitCommon/MathTools/interface/MathUtils.h"
#include "MitAna/DataUtil/interface/Debug.h"
#include "MitAna/DataTree/interface/Names.h"
#include "MitAna/DataCont/interface/ObjArray.h"
#include "MitAna/DataTree/interface/ParticleCol.h"
#include "MitPhysics/Init/interface/ModNames.h"
#include "MitPhysics/Utils/interface/JetTools.h"
#include "MitPhysics/Utils/interface/MetTools.h"

#include "MitMonoJet/SelMods/interface/BoostedVAnalysisMod.h"


using namespace mithep;
ClassImp(mithep::BoostedVAnalysisMod)

//--------------------------------------------------------------------------------------------------
  BoostedVAnalysisMod::BoostedVAnalysisMod(const char *name, const char *title) :
    BaseMod(name,title),
    // define all the Branches to load
    fMetBranchName         ("PFMet"),
    fFatJetsName           ("AKt7PFJets"),
    fJetsName              (Names::gkPFJetBrn),
    fElectronsName         (Names::gkElectronBrn),
    fMuonsName             (Names::gkMuonBrn),
    fLeptonsName           (ModNames::gkMergedLeptonsName),
    fMetFromBranch         (kTRUE),
    fFatJetsFromBranch     (kTRUE),
    fJetsFromBranch        (kTRUE),
    fElectronsFromBranch   (kTRUE),
    fMuonsFromBranch       (kTRUE),
    // define active preselection regions
    fApplyTopPresel        (kTRUE),
    fApplyWlepPresel       (kTRUE),
    fApplyZlepPresel       (kTRUE),
    fApplyMetPresel        (kTRUE), 
    fApplyVbfPresel        (kTRUE),
    // collections
    fMet                   (0),
    fJets                  (0),
    fElectrons             (0),
    fMuons                 (0),    
    // cuts
    fMinFatJetPt           (200),
    fMinTagJetPt           (100),
    fMinMet                (200),
    fMinVbfMass            (800),
    // counters
    fAll                   (0),
    fPass                  (0)
{
  // Constructor.
}

//--------------------------------------------------------------------------------------------------
void BoostedVAnalysisMod::Begin()
{
}

//--------------------------------------------------------------------------------------------------
void BoostedVAnalysisMod::SlaveBegin()
{
  // Load Branches
  ReqEventObject(fMetBranchName,   fMet,         fMetFromBranch);
  ReqEventObject(fJetsName,        fJets,        fJetsFromBranch);
  ReqEventObject(fElectronsName,   fElectrons,   fElectronsFromBranch);
  ReqEventObject(fMuonsName,       fMuons,       fMuonsFromBranch);
}

//--------------------------------------------------------------------------------------------------
void BoostedVAnalysisMod::Process()
{
  LoadEventObject(fMetBranchName,fMet,fMetFromBranch);
  LoadEventObject(fElectronsName,fElectrons,fElectronsFromBranch);
  LoadEventObject(fMuonsName,fMuons,fMuonsFromBranch);

  fJets = GetObjThisEvt<JetOArr>(fJetsName);
  fFatJets = GetObjThisEvt<JetOArr>(fFatJetsName);

  ParticleOArr   *leptons      = GetObjThisEvt<ParticleOArr>(fLeptonsName);

  // Initialize selection flags
  Bool_t passTopPresel = kFALSE;
  Bool_t passWlepPresel = kFALSE;
  Bool_t passZlepPresel = kFALSE;
  Bool_t passMetPresel = kFALSE;
  Bool_t passVbfPresel = kFALSE;

  // Increment all events counter
  fAll++;
  
  // Determine if the event passes any of the preselection cuts
  if (fApplyTopPresel) {
    // Top Preselection: require boosted jet + 2 b-jets + lepton
    // Following TOP-12-042 PAS selections
    // also require standard tag jet selection for trigger matching
    int nGoodFatJets = 0;
    int nGoodTagJets = 0;
    int nGoodBJets = 0;
    int nGoodLeptons = 0;

    // FatJets
    for (UInt_t i = 0; i < fFatJets->GetEntries(); ++i) {
      const Jet *jet = fJets->At(i);
      // Pt and eta cuts
      if (jet->Pt() < fMinFatJetPt || fabs(jet->Eta()) > 2.5)
        nGoodFatJets++;
    }

    // Jets
    for (UInt_t i = 0; i < fJets->GetEntries(); ++i) {
      const Jet *jet = fJets->At(i);
      // Pt and eta cuts
      if (jet->Pt() < 30. || fabs(jet->Eta()) > 2.5)
        continue;
      // Check high pt jet
      if (jet->Pt() > fMinTagJetPt)
        nGoodTagJets++;
      // Check b-tagging
      if (jet->CombinedSecondaryVertexBJetTagsDisc() > 0.244)
        nGoodBJets++;
    }

    // Leptons
    for (UInt_t i = 0; i < fElectrons->GetEntries(); ++i) {
      const Electron *ele = fElectrons->At(i);
      // Pt and eta cuts
      if (ele->Pt() < 30. || fabs(ele->Eta()) > 2.5 || (fabs(ele->Eta()) > 1.442 && fabs(ele->Eta()) > 1.566))
        continue;
      nGoodLeptons++;
    }
    for (UInt_t i = 0; i < fMuons->GetEntries(); ++i) {
      const Muon *mu = fMuons->At(i);
      // Pt and eta cuts
      if (mu->Pt() < 30. || fabs(mu->Eta()) > 2.1)
        continue;
      nGoodLeptons++;
    }
    if (nGoodTagJets > 0 && nGoodBJets > 1 && nGoodFatJets > 0 && nGoodLeptons > 0)
      passTopPresel = kTRUE;
  }

  if (fApplyWlepPresel) {
    // W Preselection: require high pt jet + lepton + MET
    int nGoodTagJets = 0;
    int nGoodLeptons = 0;

    // Jets
    for (UInt_t i = 0; i < fJets->GetEntries(); ++i) {
      const Jet *jet = fJets->At(i);
      // Pt and eta cuts
      if (jet->Pt() < fMinTagJetPt || fabs(jet->Eta()) > 2.5)
        continue;
      nGoodTagJets++;
    }

    // Leptons
    for (UInt_t i = 0; i < leptons->GetEntries(); ++i) {
      const Particle *lep = leptons->At(i);
      // Pt cut
      if (lep->Pt() < 10.)
        continue;
      nGoodLeptons++;
    }
    
    // Corrected MET
    float corrMetPt, corrMetPhi;
    if (leptons->GetEntries() > 0)
      CorrectMet(fMet->At(0)->Pt(), fMet->At(0)->Phi(),leptons->At(0),0,
                 corrMetPt,corrMetPhi);
        
    if (nGoodTagJets > 0 && nGoodLeptons > 0 && corrMetPt > fMinMet)
      passWlepPresel = kTRUE;
  }

  if (fApplyZlepPresel) {
    // Z Preselection: require high pt jet + di-leptons
    int nGoodTagJets = 0;
    int nGoodLeptons = 0;

    // Jets
    for (UInt_t i = 0; i < fJets->GetEntries(); ++i) {
      const Jet *jet = fJets->At(i);
      // Pt and eta cuts
      if (jet->Pt() < fMinTagJetPt || fabs(jet->Eta()) > 2.5)
        continue;
      nGoodTagJets++;
    }

    // Leptons
    for (UInt_t i = 0; i < leptons->GetEntries(); ++i) {
      const Particle *lep = leptons->At(i);
      // Pt cut
      if (lep->Pt() < 10.)
        continue;
      nGoodLeptons++;
    }
    
    // Corrected MET
    float corrMetPt, corrMetPhi;
    if (leptons->GetEntries() > 1)
      CorrectMet(fMet->At(0)->Pt(), fMet->At(0)->Phi(),leptons->At(0),leptons->At(1),
                 corrMetPt,corrMetPhi);
        
    if (nGoodTagJets > 0 && nGoodLeptons > 1 && corrMetPt > fMinMet)
      passZlepPresel = kTRUE;
  }

  if (fApplyMetPresel) {
    // Z Preselection: require boosted jet + MET
    int nGoodTagJets = 0;

    // Jets
    for (UInt_t i = 0; i < fJets->GetEntries(); ++i) {
      const Jet *jet = fJets->At(i);
      // Pt and eta cuts
      if (jet->Pt() < fMinTagJetPt || fabs(jet->Eta()) > 2.5)
        continue;
      nGoodTagJets++;
    }
        
    if (nGoodTagJets > 0 && fMet->At(0)->Pt() > fMinMet)
      passMetPresel = kTRUE;
  }

  if (fApplyVbfPresel) {
    // Vbf Preselection: require VBF jets + MET
    int nGoodVbfPairs = 0;

    // Jets
    for (UInt_t i = 0; i < fJets->GetEntries(); ++i) {
      const Jet *jetOne = fJets->At(i);
      // Pt and eta cuts
      if (jetOne->Pt() < fMinTagJetPt || fabs(jetOne->Eta()) > 2.5)
        continue;
      // Jet mulitplicity cut (at least need second jet)
      if (fJets->GetEntries() < 2) 
        continue;
      // di-jet mass cut                
      for (UInt_t j = i+1; j < fJets->GetEntries(); ++j) {
        const Jet *jetTwo = fJets->At(j);
        if ((jetOne->Mom() + jetTwo->Mom()).M() > fMinVbfMass)
          nGoodVbfPairs++;
      }
                
    }
        
    if (nGoodVbfPairs > 0 && fMet->At(0)->Pt() > fMinMet)
      passVbfPresel = kTRUE;
  }
  
  // Skip event if it does not pass any preselection
  if (!passTopPresel && !passWlepPresel && !passZlepPresel && !passMetPresel && !passVbfPresel) {
    this->SkipEvent(); 
    return;
  }
   
  // Increment passed events counter
  fPass++;
  return;
}

//--------------------------------------------------------------------------------------------------
void BoostedVAnalysisMod::SlaveTerminate()
{

  // Terminate preselection and print out module selection efficiency
  Double_t frac =  100.*fPass/fAll;
  Info("SlaveTerminate", "Selected %.2f%% events (%lld out of %lld)", 
       frac, fPass, fAll);
  
  return;
}

//--------------------------------------------------------------------------------------------------
void BoostedVAnalysisMod::Terminate()
{
}

//--------------------------------------------------------------------------------------------------
void BoostedVAnalysisMod::CorrectMet(const float met, const float metPhi,
                                     const Particle *l1, const Particle *l2,
                                     float &newMet, float &newMetPhi)
{
  // inputs: met, metPhi, l1, [ l2 only used if pointer is non-zero ]
  // outputs: newMet, newMetPhi

  float newMetX;
  float newMetY;
  if (l2) { // these are doubles on the right: need full calculation in one line
    newMetX = met*TMath::Cos(metPhi) + l1->Mom().Px() + l2->Mom().Px();
    newMetY = met*TMath::Sin(metPhi) + l1->Mom().Py() + l2->Mom().Py();
  }
  else {
    newMetX = met*TMath::Cos(metPhi) + l1->Mom().Px();
    newMetY = met*TMath::Sin(metPhi) + l1->Mom().Py();
  }
  
  newMet = TMath::Sqrt(TMath::Power(newMetX,2) + TMath::Power(newMetY,2));
  newMetPhi = TMath::ATan2(newMetY,newMetX);
}


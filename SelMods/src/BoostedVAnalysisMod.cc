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
    fPhotonsName           (Names::gkPhotonBrn),
    fLeptonsName           (ModNames::gkMergedLeptonsName),
    fMetFromBranch         (kTRUE),
    fFatJetsFromBranch     (kTRUE),
    fJetsFromBranch        (kTRUE),
    fElectronsFromBranch   (kTRUE),
    fMuonsFromBranch       (kTRUE),
    fPhotonsFromBranch     (kTRUE),
    // define active preselection regions
    fApplyTopPresel        (kTRUE),
    fApplyWlepPresel       (kTRUE),
    fApplyZlepPresel       (kTRUE),
    fApplyMetPresel        (kTRUE), 
    fApplyVbfPresel        (kTRUE),
    fApplyGjetPresel       (kTRUE),
    fApplyFatJetPresel     (kTRUE),
    fFillAndPublishPresel  (kTRUE),
    // collections
    fMet                   (0),
    fJets                  (0),
    fElectrons             (0),
    fMuons                 (0),
    fPhotons               (0),
    fEvtSelData            (0),    
    // cuts
    fMinFatJetPt           (200),
    fMinTagJetPt           (100),
    fMinVbfJetPt           (40),
    fMinMet                (200),
    fMinVbfMass            (800),
    fMinVbfMet             (110),
    fMinPhotonPt           (150),
    // counters
    fAll                   (0),
    fPass                  (0)
{
  // Constructor.
}

//--------------------------------------------------------------------------------------------------
BoostedVAnalysisMod::~BoostedVAnalysisMod()
{
  // Destructor
  if (fXlEvtSelData)
    delete fXlEvtSelData;
}

//--------------------------------------------------------------------------------------------------
void BoostedVAnalysisMod::Begin()
{
}

//--------------------------------------------------------------------------------------------------
void BoostedVAnalysisMod::SlaveBegin()
{
  // Load Branches
  ReqEventObject(fMetBranchName, fMet,       fMetFromBranch);
  ReqEventObject(fJetsName,      fJets,      fJetsFromBranch);
  ReqEventObject(fElectronsName, fElectrons, fElectronsFromBranch);
  ReqEventObject(fMuonsName,     fMuons,     fMuonsFromBranch);
  ReqEventObject(fPhotonsName,   fPhotons,   fPhotonsFromBranch);
  
  // If requested by the user prepare the collections to store the
  // preselection word
  if (fFillAndPublishPresel) {
    ReqEventObject(Names::gkEvtSelDataBrn, fEvtSelData, true);
    // Create the new output collection
    fXlEvtSelData = new XlEvtSelData("XlEvtSelData");
    PublishObj(fXlEvtSelData);
  }

}

//--------------------------------------------------------------------------------------------------
void BoostedVAnalysisMod::Process()
{
  LoadEventObject(fMetBranchName, fMet,       fMetFromBranch);
  LoadEventObject(fElectronsName, fElectrons, fElectronsFromBranch);
  LoadEventObject(fMuonsName,     fMuons,     fMuonsFromBranch);
  LoadEventObject(fPhotonsName,   fPhotons,   fPhotonsFromBranch);
                 
  if (fFillAndPublishPresel) {
    LoadEventObject(Names::gkEvtSelDataBrn, fEvtSelData, true);
  }

  fJets    = GetObjThisEvt<JetOArr>(fJetsName);
  if (fApplyFatJetPresel)
    fFatJets = GetObjThisEvt<JetOArr>(fFatJetsName);

  ParticleOArr *leptons = GetObjThisEvt<ParticleOArr>(fLeptonsName);
  
  // Initialize selection flags
  Bool_t passTopPresel = kFALSE;
  Bool_t passWlepPresel = kFALSE;
  Bool_t passZlepPresel = kFALSE;
  Bool_t passMetPresel = kFALSE;
  Bool_t passVbfPresel = kFALSE;
  Bool_t passGjetPresel = kFALSE;

  // Increment all events counter
  fAll++;
  
  // Discard events with no jets
  if (fJets->GetEntries() < 1) {
    this->SkipEvent(); 
    return;
  }

  // Determine if the event passes any of the preselection cuts
  if (fApplyTopPresel) {
    // Top Preselection: require boosted jet + 2 b-jets + lepton
    // Following TOP-12-042 PAS selections
    // also require standard tag jet selection for trigger matching
    int nGoodFatJets = 0;
    int nGoodTagJets = 0;
    int nGoodBJets = 0;
    int nGoodLeptons = 0;

    // FatJets: if fatJets not available adjust preselection params
    if (!fApplyFatJetPresel)
      nGoodFatJets = 1;
    if (fApplyFatJetPresel && fFatJets->GetEntries() > 0) { 
      for (UInt_t i = 0; i < fFatJets->GetEntries(); ++i) {
        const Jet *jet = fFatJets->At(i);
        // Pt and eta cuts
        if (jet->Pt() < fMinFatJetPt || fabs(jet->Eta()) > 2.5)
          nGoodFatJets++;
      }
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
    if (fElectrons->GetEntries() > 0) {
      for (UInt_t i = 0; i < fElectrons->GetEntries(); ++i) {
        const Electron *ele = fElectrons->At(i);
        // Pt and eta cuts
        if (ele->Pt() < 30. || fabs(ele->Eta()) > 2.5 || (fabs(ele->Eta()) > 1.442 && fabs(ele->Eta()) > 1.566))
          continue;
        nGoodLeptons++;
      }
    }
    if (fMuons->GetEntries() > 0) {
      for (UInt_t i = 0; i < fMuons->GetEntries(); ++i) {
        const Muon *mu = fMuons->At(i);
        // Pt and eta cuts
        if (mu->Pt() < 30. || fabs(mu->Eta()) > 2.1)
          continue;
        nGoodLeptons++;
      }
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
    if (leptons->GetEntries() > 0) {
      for (UInt_t i = 0; i < leptons->GetEntries(); ++i) {
        const Particle *lep = leptons->At(i);
        // Pt cut
        if (lep->Pt() < 10.)
          continue;
        nGoodLeptons++;
      }
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
    if (leptons->GetEntries() > 0) {
      for (UInt_t i = 0; i < leptons->GetEntries(); ++i) {
        const Particle *lep = leptons->At(i);
        // Pt cut
        if (lep->Pt() < 10.)
          continue;
        nGoodLeptons++;
      }
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
    // Vbf Preselection: require VBF jets + METnoMu
    int nGoodVbfPairs = 0;

    // Jets
    for (UInt_t i = 0; i < fJets->GetEntries(); ++i) {
      const Jet *jetOne = fJets->At(i);
      // Pt and eta cuts
      if (jetOne->Pt() < fMinVbfJetPt)
        continue;
      // Jet mulitplicity cut (at least need second jet)
      if (fJets->GetEntries() < 2) 
        continue;
      // di-jet mass cut                
      for (UInt_t j = i+1; j < fJets->GetEntries(); ++j) {
        const Jet *jetTwo = fJets->At(j);
        if (jetTwo->Pt() < fMinVbfJetPt) 
          continue;
        if ((jetOne->Mom() + jetTwo->Mom()).M() > fMinVbfMass)
          nGoodVbfPairs++;
      }
                
    }
    
    // METnoMu
    float METnoMuPt, METnoMuPhi;
    if (fMuons->GetEntries() > 1) {
      CorrectMet(fMet->At(0)->Pt(), fMet->At(0)->Phi(),fMuons->At(0),fMuons->At(1),
                 METnoMuPt,METnoMuPhi);
    }
    else if (fMuons->GetEntries() == 1) {
      CorrectMet(fMet->At(0)->Pt(), fMet->At(0)->Phi(),fMuons->At(0),0,
                 METnoMuPt,METnoMuPhi);
    }
    else 
      METnoMuPt = fMet->At(0)->Pt();                      
        
    if (nGoodVbfPairs > 0 && METnoMuPt > fMinVbfMet)
      passVbfPresel = kTRUE;
  }

  if (fApplyGjetPresel) {
    // G+jets Preselection: require boosted jet + photon
    int nGoodFatJets = 0;
    int nGoodPhotons = 0;

    // FatJets: if fatJets not available adjust preselection params
    if (!fApplyFatJetPresel)
      nGoodFatJets = 1;
    if (fApplyFatJetPresel && fFatJets->GetEntries() > 0) { 
      for (UInt_t i = 0; i < fFatJets->GetEntries(); ++i) {
        const Jet *jet = fFatJets->At(i);
        // Pt and eta cuts
        if (jet->Pt() < fMinFatJetPt || fabs(jet->Eta()) > 2.5)
          nGoodFatJets++;
      }
    }

    // Photons
    if (fPhotons->GetEntries() > 0) {
      for (UInt_t i = 0; i < fPhotons->GetEntries(); ++i) {
        const Photon *pho = fPhotons->At(i);
        // Pt and eta cuts
        if (pho->Pt() < fMinPhotonPt || fabs(pho->Eta()) > 2.5)
          continue;
        nGoodPhotons++;
      }
    }
    if (nGoodFatJets > 0 && nGoodPhotons > 0)
      passGjetPresel = kTRUE;
  }

  // Skip event if it does not pass any preselection
  if (!passTopPresel && !passWlepPresel && !passZlepPresel 
   && !passMetPresel && !passVbfPresel  && !passGjetPresel) {
    this->SkipEvent(); 
    return;
  }

  // Store the preselection word in the extended event selection object
  if (fFillAndPublishPresel) {
    int preselectionWord = GetPreselWord (
                           passTopPresel, 
                           passWlepPresel,
                           passZlepPresel,
                           passMetPresel, 
                           passVbfPresel, 
                           passGjetPresel); 
    fXlEvtSelData->SetFiltersWord(fEvtSelData->metFiltersWord());
    fXlEvtSelData->SetPreselWord(preselectionWord);
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

//--------------------------------------------------------------------------------------------------
int  BoostedVAnalysisMod::GetPreselWord( 
                          bool passTopPresel, 
                          bool passWlepPresel,
                          bool passZlepPresel,
                          bool passMetPresel, 
                          bool passVbfPresel,
                          bool passGjetPresel)
{  
  // This function creates the word containing the bit decisions.
  // The bit ordering follows the order of the parameters passed
  // to this function. 
  
  //Initialize the word
  int theWord = 0;
  //Initialize the vector of bits
  std::vector<int> theBits;
  theBits.push_back((int) passTopPresel);
  theBits.push_back((int) passWlepPresel);
  theBits.push_back((int) passZlepPresel);
  theBits.push_back((int) passMetPresel);
  theBits.push_back((int) passVbfPresel);
  theBits.push_back((int) passGjetPresel);
  //Create the word
  for (unsigned int iBit = 0; iBit < theBits.size(); iBit++)
    theWord |= theBits[iBit] << iBit;
  
  return theWord;
}
  

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
#include "MitAna/DataTree/interface/MCParticleFwd.h"
#include "MitPhysics/Init/interface/ModNames.h"
#include "MitPhysics/Utils/interface/JetTools.h"
#include "MitPhysics/Utils/interface/MetTools.h"

#include "MitMonoJet/SelMods/interface/MonoJetAnalysisMod.h"


using namespace mithep;
ClassImp(mithep::MonoJetAnalysisMod)

//--------------------------------------------------------------------------------------------------
  MonoJetAnalysisMod::MonoJetAnalysisMod(const char *name, const char *title) :
    BaseMod(name,title),
    //define all the Branches to load
    fCategoriesName        ("MonoJetCategories"),
    fMetBranchName         ("PFMet"),
    fJetsName              (Names::gkPFJetBrn),
    fElectronsName         (Names::gkElectronBrn),
    fMuonsName             (Names::gkMuonBrn),
    fTausName              (Names::gkPFTauBrn),
    fLeptonsName           (ModNames::gkMergedLeptonsName),
    fPFCandidatesName      (Names::gkPFCandidatesBrn),
    fMetFromBranch         (kTRUE),
    fJetsFromBranch        (kTRUE),
    fPFCandidatesFromBranch(kTRUE),
    fElectronsFromBranch   (kTRUE),
    fMuonsFromBranch       (kTRUE),
    fTausFromBranch        (kTRUE),
    // collections
    fMet                   (0),
    fJets                  (0),
    fElectrons             (0),
    fMuons                 (0),
    fPFTaus                (0),
    fPFCandidates          (0),
    // event category flags
    fCategories            (0),
    // counters
    fNEventsSelected       (0)
{
  // cuts
  std::fill_n(fMinNumLeptons, nCat, 0);
  std::fill_n(fMinNumTaus, nCat, 0);
  std::fill_n(fMinNumJets, nCat, 1);
  std::fill_n(fMinNumGenNeutrinos, nCat, 0);
  std::fill_n(fMinJetEt, nCat, 30.);
  std::fill_n(fMaxJetEta, nCat, 2.4);
  std::fill_n(fMinMetEt, nCat, 30.);
  std::fill_n(fMinChargedHadronFrac, nCat, 0.);
  std::fill_n(fMaxNeutralHadronFrac, nCat, 1.);
  std::fill_n(fMaxNeutralEmFrac, nCat, 1.);
}

//--------------------------------------------------------------------------------------------------
void MonoJetAnalysisMod::Begin()
{
}

//--------------------------------------------------------------------------------------------------
void MonoJetAnalysisMod::SlaveBegin()
{
  // Load Branches
  ReqEventObject(fMetBranchName,   fMet,         fMetFromBranch);
  ReqEventObject(fJetsName,        fJets,        fJetsFromBranch);
  ReqEventObject(fPFCandidatesName,fPFCandidates,fPFCandidatesFromBranch);
  ReqEventObject(fElectronsName,   fElectrons,   fElectronsFromBranch);
  ReqEventObject(fMuonsName,       fMuons,       fMuonsFromBranch);
  ReqEventObject(fTausName,        fPFTaus,      fTausFromBranch);

  // Selection Histograms
  for(unsigned iCat = 0; iCat != nCat; ++iCat){
    AddTH1(fMonoJetSelection[iCat],TString::Format("hMonoJetSelection%d", iCat), ";Cut Number;Number of Events",17,-1.5,15.5);
    // Create const char* labels with cut values
    std::ostringstream os1;
    std::string s1;
    os1<<"N Jets >= "<<fMinNumJets[iCat];
    s1 = os1.str();
    const char* labelCut1 = s1.c_str();
    std::ostringstream os2;
    std::string s2;
    os2<<"Jet Et >= "<<fMinJetEt[iCat];
    s2 = os2.str();
    const char* labelCut2 = s2.c_str();
    std::ostringstream os3;
    std::string s3;
    os3<<"Jet Eta <= "<<fMaxJetEta[iCat];
    s3 = os3.str();
    const char* labelCut3 = s3.c_str();
    std::ostringstream os4;
    std::string s4;
    os4<<"Met >= "<<fMinMetEt[iCat];
    s4 = os4.str();
    const char* labelCut4 = s4.c_str();
    std::ostringstream os5;
    std::string s5;
    os5<<"Charged Hadron Fraction >= "<<fMinChargedHadronFrac[iCat];
    s5 = os5.str();
    const char* labelCut5 = s5.c_str();
    std::ostringstream os6;
    std::string s6;
    os6<<"Neutral Hadron Fraction <= "<<fMaxNeutralHadronFrac[iCat];
    s6 = os6.str();
    const char* labelCut6 = s6.c_str();
    std::ostringstream os7;
    std::string s7;
    os7<<"Neutral Em Fraction <= "<<fMaxNeutralEmFrac[iCat];
    s7 = os7.str();
    const char* labelCut7 = s7.c_str();
    std::ostringstream os8;
    std::string s8;
    os8<<"N Leptons >= "<<fMinNumLeptons[iCat];
    s8 = os8.str();
    const char* labelCut8 = s8.c_str();
    std::ostringstream os9;
    std::string s9;
    os9<<"N Gen Neutrinos >= "<<fMinNumGenNeutrinos[iCat];
    s9 = os9.str();
    const char* labelCut9 = s9.c_str();

    // Set selection histogram bin labels
    fMonoJetSelection[iCat]->GetXaxis()->TAxis::SetBinLabel(1, "All Events");
    fMonoJetSelection[iCat]->GetXaxis()->TAxis::SetBinLabel(2, labelCut1);
    fMonoJetSelection[iCat]->GetXaxis()->TAxis::SetBinLabel(3, labelCut2);
    fMonoJetSelection[iCat]->GetXaxis()->TAxis::SetBinLabel(4, labelCut3);
    fMonoJetSelection[iCat]->GetXaxis()->TAxis::SetBinLabel(5, labelCut4);
    fMonoJetSelection[iCat]->GetXaxis()->TAxis::SetBinLabel(6, labelCut5);
    fMonoJetSelection[iCat]->GetXaxis()->TAxis::SetBinLabel(7, labelCut6);
    fMonoJetSelection[iCat]->GetXaxis()->TAxis::SetBinLabel(8, labelCut7);
    fMonoJetSelection[iCat]->GetXaxis()->TAxis::SetBinLabel(9, labelCut8);
    fMonoJetSelection[iCat]->GetXaxis()->TAxis::SetBinLabel(10, labelCut9);

    // Histograms after preselection
    AddTH1(fMetEt[iCat], TString::Format("hMetEt%d", iCat), ";MetEt;Number of Events" ,400,0.,400.);
    AddTH1(fJetEt[iCat], TString::Format("hJetEt%d", iCat), ";JetEt;Number of Events" ,400,0.,400.);
    AddTH1(fJetEta[iCat],TString::Format("hJetEta%d", iCat),";JetEta;Number of Events",50,-5,5);
  }

  fCategories = new mithep::Array<mithep::TriggerMask>(1, fCategoriesName);
  fCategories->AddNew();
  fCategories->SetMustDeleteBit();
  PublishObj(fCategories);
}

//--------------------------------------------------------------------------------------------------
void MonoJetAnalysisMod::Process()
{
  LoadEventObject(fMetBranchName,fMet,fMetFromBranch);
  LoadEventObject(fElectronsName,fElectrons,fElectronsFromBranch);
  LoadEventObject(fMuonsName,fMuons,fMuonsFromBranch);
  LoadEventObject(fTausName,fPFTaus,fTausFromBranch);
  LoadEventObject(fPFCandidatesName,fPFCandidates,fPFCandidatesFromBranch);

  fJets = GetObjThisEvt<JetOArr>(fJetsName);

  ParticleOArr   *leptons      = GetObjThisEvt<ParticleOArr>(fLeptonsName);
  MCParticleOArr *genNeutrinos = GetObjThisEvt<MCParticleOArr>(ModNames::gkMCNeutrinosName);

  // Define Cuts
  const int nCuts = 10;
  bool passCut[nCat][nCuts];
  for(unsigned iCat = 0; iCat != nCat; ++iCat)
    std::fill_n(passCut[iCat], nCuts, false);

  // Discard events with no identified required objects (jets, leptons, taus)
  for(unsigned iCat = 0; iCat != nCat; ++iCat){
    if (fJets->GetEntries() >= fMinNumJets[iCat])
      passCut[iCat][0] = true;
    if (leptons->GetEntries() >= fMinNumLeptons[iCat])
      passCut[iCat][1] = true;
    if (fPFTaus->GetEntries() >= fMinNumTaus[iCat])
      passCut[iCat][2] = true;
  }

  // Discard events with soft jets
  for(unsigned iCat = 0; iCat != nCat; ++iCat){
    int nHardJet = 0;
    for (UInt_t i = 0; i < fJets->GetEntries(); ++i) {
      const Jet *jet = fJets->At(i);
      if (jet->Et() <= fMinJetEt[iCat])
        continue;
      nHardJet++;
    }
    if (nHardJet > 0)
      passCut[iCat][3] = true;
  }

  // Discard events that are outside of eta range
  vector<int> theGoodJets[nCat];
  for(unsigned iCat = 0; iCat != nCat; ++iCat){
    int nHardEtaJet = 0;
    for (UInt_t i = 0; i < fJets->GetEntries(); ++i) {
      const Jet *jet = fJets->At(i);
      if (jet->Et() <= fMinJetEt[iCat])
        continue;
      if (TMath::Abs(jet->Eta()) >= fMaxJetEta[iCat])
        continue;
      nHardEtaJet++;
      theGoodJets[iCat].push_back((int) i);
    }
    if (nHardEtaJet > 0){
      passCut[iCat][4] = true;
    }
  }

  // Discard events with soft met
  const Met *stdMet = fMet->At(0);
  for(unsigned iCat = 0; iCat != nCat; ++iCat){
    if ( stdMet->Pt() > fMinMetEt[iCat] )
    passCut[iCat][5] = true;
  }

  for(unsigned iCat = 0; iCat != nCat; ++iCat){
    passCut[iCat][6]=true;
    passCut[iCat][7]=true;
    passCut[iCat][8]=true;
  }

  // Make Selection Histograms. Number of events passing each level of cut

  // Cut Selection Histograms

  for(unsigned iCat = 0; iCat != nCat; ++iCat){
    fMonoJetSelection[iCat]->Fill(-1,1);

    for (int k=0;k<nCuts;k++) {
      bool pass = true;
      bool passPreviousCut = true;
      for (int p=0;p<=k;p++) {
        pass = (pass && passCut[iCat][p]);
        if (p<k)
          passPreviousCut = (passPreviousCut&& passCut[iCat][p]);
      }
      if (pass)
        fMonoJetSelection[iCat]->Fill(k,1);
    }
  }

  UInt_t nuInAcceptance = 0;
  for (UInt_t i = 0; i < genNeutrinos->GetEntries(); ++i)
    if (fabs(genNeutrinos->At(i)->Eta())<2.1 and genNeutrinos->At(i)->Pt()>10) nuInAcceptance++;

  for(unsigned iCat = 0; iCat != nCat; ++iCat){
    if(nuInAcceptance >= fMinNumGenNeutrinos[iCat])
      passCut[iCat][9] = true;
  }

  mithep::BitMask1024 catBits;
  for(unsigned iCat = 0; iCat != nCat; ++iCat){
    bool passAllCuts = true;
    for (int c=0; c<nCuts; c++) passAllCuts = passAllCuts & passCut[iCat][c];
    if (passAllCuts) {
      fNEventsSelected++;
      // Make Preselection Histograms
      for (int iGoodJet=0; iGoodJet<(int)theGoodJets[iCat].size(); iGoodJet++)
        fJetEt[iCat]->Fill(fJets->At(theGoodJets[iCat][iGoodJet])->Et());
      fMetEt[iCat]->Fill(fMet->At(0)->Et());
      fJetEta[iCat]->Fill(fJets->At(0)->Eta());

      catBits.SetBit(iCat);
    }
  }

  if(catBits.NBitsSet() == 0){
    this->SkipEvent(); // skip the event if does not passes the cuts
  }

  fCategories->At(0)->SetBits(catBits);
}

//--------------------------------------------------------------------------------------------------
void MonoJetAnalysisMod::SlaveTerminate()
{
  cout << "selected events on MonoJetAnalysisMod: " << fNEventsSelected << endl;

  fCategories->Delete();
}
//--------------------------------------------------------------------------------------------------
void MonoJetAnalysisMod::Terminate()
{
}

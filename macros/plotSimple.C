#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"

void plotSimple(){
  TFile* inFile = new TFile("boostedv_r12b-smu-j22-v1_noskim_000.root");
  TTree* inTree = (TTree*) inFile->Get("Events");
  
  TCanvas *can = new TCanvas();
  inTree->Draw("XlFatJets->At(0)->SubJet(0)->Pt()","XlSubJets->GetEntries()>0");

  TCanvas *can2 = new TCanvas();
  inTree->Draw("XlFatJets->At(0)->SubJet(0)->SubJetType()","XlSubJets->GetEntries()>0");

  TCanvas *can3 = new TCanvas();
  inTree->Draw("XlFatJets->At(0)->SubJet(2)->SubJetType()","XlSubJets->GetEntries()>2");

  TCanvas *can4 = new TCanvas();
  inTree->Draw("XlFatJets->At(0)->Tau1()","XlSubJets->GetEntries()>2");

  TCanvas *can5 = new TCanvas();
  inTree->Draw("XlFatJets->At(0)->ChargedHadronEnergy()/XlFatJets->At(0)->RawMom().E();","XlSubJets->GetEntries()>2");

  TCanvas *can6 = new TCanvas();
  inTree->Draw("XlFatJets->At(0)->GroomedMom().M()","XlSubJets->GetEntries()>2");
  
  return;
}

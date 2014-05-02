#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"

void plotSimple(){
  TFile* inFile = new TFile("skimtest_000.root");
  TTree* inTree = (TTree*) inFile->Get("Events");
  
  TCanvas *can = new TCanvas();
  inTree->Draw("XlFatJets->At(0)->SubJet(0)->Pt()","XlSubJets->GetEntries()>0");

  TCanvas *can2 = new TCanvas();
  inTree->Draw("XlFatJets->At(0)->SubJet(0)->SubJetType()","XlSubJets->GetEntries()>0");

  TCanvas *can3 = new TCanvas();
  inTree->Draw("XlFatJets->At(0)->SubJet(2)->SubJetType()","XlSubJets->GetEntries()>2");

  TCanvas *can4 = new TCanvas();
  inTree->Draw("XlFatJets->At(0)->Tau1()","XlSubJets->GetEntries()>2");
  
  return;
}

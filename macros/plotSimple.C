#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"

void plotSimple(){
  //TFile* inFile = new TFile("/scratch4/dimatteo/cms/hist/boostedv-v3/t2mit/filefi/032/s12-pj170_300-v7a/boostedv-v3_s12-pj170_300-v7a_noskim_0006_ntuple_000.root");
  TFile* inFile = new TFile("/scratch1/dimatteo/cmssw/033/CMSSW_5_3_16/src/MitMonoJet/macros/MonoVSynch/boostedv_s12-ttj-sl-v1-v7c_noskim_0000_ntuple_000.root");
  TTree* inTree = (TTree*) inFile->Get("Events");

  TCanvas *can = new TCanvas();
  inTree->Draw("XsMuons->At(0)->ChHadIso()","XsMuons->GetEntries()>0");

  TCanvas *can2 = new TCanvas();
  inTree->Draw("SkmCleanJets->At(0)->NPFCands()","XlSubJets->GetEntries()>0");

  TCanvas *can3 = new TCanvas();
  inTree->Draw("XsTaus->At(0)->Mom().Pt()","XsTaus->GetEntries()>0");
  
  //TCanvas *can = new TCanvas();
  //inTree->Draw("(XlFatJets->At(0)->Mom().Pt()-SkmCleanJets->At(0)->Mom().Pt())/SkmCleanJets->At(0)->Mom().Pt()","XlSubJets->GetEntries()>0");

  //TCanvas *can2 = new TCanvas();
  //inTree->Draw("XlFatJets->At(0)->SubJet(0)->SubJetType()","XlSubJets->GetEntries()>0");

  //TCanvas *can3 = new TCanvas();
  //inTree->Draw("XlFatJets->At(0)->SubJet(2)->SubJetType()","XlSubJets->GetEntries()>2");

  //TCanvas *can4 = new TCanvas();
  //inTree->Draw("XlFatJets->At(0)->Tau1()","XlSubJets->GetEntries()>2");

  //TCanvas *can5 = new TCanvas();
  //inTree->Draw("XlFatJets->At(0)->ChargedHadronEnergy()/XlFatJets->At(0)->RawMom().E();","XlSubJets->GetEntries()>2");

  //TCanvas *can7 = new TCanvas();
  //inTree->Draw("XlFatJets->At(0)->QJetVol()","XlSubJets->GetEntries()>2");

  //TCanvas *can8 = new TCanvas();
  //inTree->Draw("XlFatJets->At(0)->C2b0()","XlSubJets->GetEntries()>2");

  //TCanvas *can9 = new TCanvas();
  //inTree->Draw("XlFatJets->At(0)->C2b0p2()","XlSubJets->GetEntries()>2");

  //TCanvas *can10 = new TCanvas();
  //inTree->Draw("XlFatJets->At(0)->C2b0p5()","XlSubJets->GetEntries()>2");

  //TCanvas *can11 = new TCanvas();
  //inTree->Draw("XlFatJets->At(0)->C2b1()","XlSubJets->GetEntries()>2");

  //TCanvas *can12 = new TCanvas();
  //inTree->Draw("XlFatJets->At(0)->C2b2()","XlSubJets->GetEntries()>2");

  //TCanvas *can13 = new TCanvas();
  //inTree->Draw("XlFatJets->At(0)->MassSDb0()","XlSubJets->GetEntries()>2");

  //TCanvas *can14 = new TCanvas();
  //inTree->Draw("XlFatJets->At(0)->MassSDb2()","XlSubJets->GetEntries()>2");

  //TCanvas *can15 = new TCanvas();
  //inTree->Draw("XlFatJets->At(0)->MassSDbm1()","XlSubJets->GetEntries()>2");

  //TCanvas *can16 = new TCanvas();
  //inTree->Draw("XlFatJets->At(0)->MassPruned()","XlSubJets->GetEntries()>2");

  //TCanvas *can17 = new TCanvas();
  //inTree->Draw("XlFatJets->At(0)->MassFiltered()","XlSubJets->GetEntries()>2");

  //TCanvas *can18 = new TCanvas();
  //inTree->Draw("XlFatJets->At(0)->MassTrimmed()","XlSubJets->GetEntries()>2");

  //TCanvas *can19 = new TCanvas();
  //inTree->Draw("XlFatJets->At(0)->QJetVol()","XlSubJets->GetEntries()>2");

  //TCanvas *can20 = new TCanvas();
  ////inTree->Draw("mithep::XlEvtSelData->metFiltersWord()","XlSubJets->GetEntries()>2");
  //inTree->Draw("mithep::XlEvtSelData->preselWord()","XlSubJets->GetEntries()>2");
  
  return;
}

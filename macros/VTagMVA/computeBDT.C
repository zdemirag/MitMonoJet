#include <cstdlib>
#include <iostream>
#include <map>
#include <string>

#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TTreeFormula.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TROOT.h"

#include "TMVA/TMVAGui.C"

#if not defined(__CINT__) || defined(__MAKECINT__)
// needs to be included when makecint runs (ACLIC)
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCategory.h"
#include "TMVA/Tools.h"
#endif

void computeBDT(std::string iName="../samples/boostedv-v8_s12-zll-ptz100-v7c_noskim_flatntuple.root",
              std::string iWeightFile="weights/TMVAClassificationCategory_BDT_simple_alpha.weights.xml") { 
  TMVA::Tools::Instance();
  TMVA::Reader *reader = new TMVA::Reader( "!Color:!Silent" );    

  float jet1QGtagComb = 0.; reader->AddVariable("2*fjet1QGtagSub2+fjet1QGtagSub1",&jet1QGtagComb);
  float fjet1QGtagSub1 = 0.; reader->AddVariable("fjet1QGtagSub1",&fjet1QGtagSub1);
  float fjet1QGtagSub2 = 0.; reader->AddVariable("fjet1QGtagSub2",&fjet1QGtagSub2);
  float fjet1QGtag = 0.; reader->AddVariable("fjet1QGtag",&fjet1QGtag);
  float fjet1PullAngle = 0.; reader->AddVariable("fjet1PullAngle",&fjet1PullAngle);
  float fjet1Pull = 0.; reader->AddVariable("fjet1Pull",&fjet1Pull);
  float fjet1MassTrimmed = 0.; reader->AddVariable("fjet1MassTrimmed",&fjet1MassTrimmed);
  float fjet1MassPruned = 0.; reader->AddVariable("fjet1MassPruned",&fjet1MassPruned);
  float fjet1MassSDbm1 = 0.; reader->AddVariable("fjet1MassSDbm1",&fjet1MassSDbm1);
  float fjet1MassSDb2 = 0.; reader->AddVariable("fjet1MassSDb2",&fjet1MassSDb2);
  float fjet1MassSDb0 = 0.; reader->AddVariable("fjet1MassSDb0",&fjet1MassSDb0);
  float fjet1QJetVol = 0.; reader->AddVariable("fjet1QJetVol",&fjet1QJetVol);
  float fjet1C2b2 = 0.; reader->AddVariable("fjet1C2b2",&fjet1C2b2);
  float fjet1C2b1 = 0.; reader->AddVariable("fjet1C2b1",&fjet1C2b1);
  float fjet1C2b0p5 = 0.; reader->AddVariable("fjet1C2b0p5",&fjet1C2b0p5);
  float fjet1C2b0p2 = 0.; reader->AddVariable("fjet1C2b0p2",&fjet1C2b0p2);
  float fjet1C2b0 = 0.; reader->AddVariable("fjet1C2b0",&fjet1C2b0);
  float fjet1Tau2 = 0.; reader->AddVariable("fjet1Tau2",&fjet1Tau2);
  float fjet1Tau1 = 0.; reader->AddVariable("fjet1Tau1",&fjet1Tau1);
  float tau2tau1 = 0.; reader->AddVariable("fjet1Tau2/fjet1Tau1",&tau2tau1);

  std::string lJetName = "BDT";
  reader->BookMVA(lJetName .c_str(),iWeightFile.c_str());
  
  TFile *lFile = new TFile(iName.c_str());
  TTree *lTree = (TTree*) lFile->FindObjectAny("DMSTree");
  TString pExpress0 = "2*fjet1QGtagSub2+fjet1QGtagSub1";
  TString pExpr0(pExpress0);
  TTreeFormula* lFVars0 = new TTreeFormula(pExpr0,pExpr0,lTree);
  TString pExpress1 = "fjet1QGtagSub1";
  TString pExpr1(pExpress1);
  TTreeFormula* lFVars1 = new TTreeFormula(pExpr1,pExpr1,lTree);
  TString pExpress2 = "fjet1QGtagSub2";
  TString pExpr2(pExpress2);
  TTreeFormula* lFVars2 = new TTreeFormula(pExpr2,pExpr2,lTree);
  TString pExpress3 = "fjet1QGtag";
  TString pExpr3(pExpress3);
  TTreeFormula* lFVars3 = new TTreeFormula(pExpr3,pExpr3,lTree);
  TString pExpress4 = "fjet1PullAngle";
  TString pExpr4(pExpress4);
  TTreeFormula* lFVars4 = new TTreeFormula(pExpr4,pExpr4,lTree);
  TString pExpress5 = "fjet1Pull";
  TString pExpr5(pExpress5);
  TTreeFormula* lFVars5 = new TTreeFormula(pExpr5,pExpr5,lTree);
  TString pExpress6 = "fjet1MassTrimmed";
  TString pExpr6(pExpress6);
  TTreeFormula* lFVars6 = new TTreeFormula(pExpr6,pExpr6,lTree);
  TString pExpress7 = "fjet1MassPruned";
  TString pExpr7(pExpress7);
  TTreeFormula* lFVars7 = new TTreeFormula(pExpr7,pExpr7,lTree);
  TString pExpress8 = "fjet1MassSDbm1";
  TString pExpr8(pExpress8);
  TTreeFormula* lFVars8 = new TTreeFormula(pExpr8,pExpr8,lTree);
  TString pExpress9 = "fjet1MassSDb2";
  TString pExpr9(pExpress9);
  TTreeFormula* lFVars9 = new TTreeFormula(pExpr9,pExpr9,lTree);
  TString pExpress10 = "fjet1MassSDb0";
  TString pExpr10(pExpress10);
  TTreeFormula* lFVars10 = new TTreeFormula(pExpr10,pExpr10,lTree);
  TString pExpress11 = "fjet1QJetVol";
  TString pExpr11(pExpress11);
  TTreeFormula* lFVars11 = new TTreeFormula(pExpr11,pExpr11,lTree);
  TString pExpress12 = "fjet1C2b2";
  TString pExpr12(pExpress12);
  TTreeFormula* lFVars12 = new TTreeFormula(pExpr12,pExpr12,lTree);
  TString pExpress13 = "fjet1C2b1";
  TString pExpr13(pExpress13);
  TTreeFormula* lFVars13 = new TTreeFormula(pExpr13,pExpr13,lTree);
  TString pExpress14 = "fjet1C2b0p5";
  TString pExpr14(pExpress14);
  TTreeFormula* lFVars14 = new TTreeFormula(pExpr14,pExpr14,lTree);
  TString pExpress15 = "fjet1C2b0p2";
  TString pExpr15(pExpress15);
  TTreeFormula* lFVars15 = new TTreeFormula(pExpr15,pExpr15,lTree);
  TString pExpress16 = "fjet1C2b0";
  TString pExpr16(pExpress16);
  TTreeFormula* lFVars16 = new TTreeFormula(pExpr16,pExpr16,lTree);
  TString pExpress17 = "fjet1Tau2";
  TString pExpr17(pExpress17);
  TTreeFormula* lFVars17 = new TTreeFormula(pExpr17,pExpr17,lTree);
  TString pExpress18 = "fjet1Tau1";
  TString pExpr18(pExpress18);
  TTreeFormula* lFVars18 = new TTreeFormula(pExpr18,pExpr18,lTree);
  TString pExpress19 = "fjet1Tau2/fjet1Tau1";
  TString pExpr19(pExpress19);
  TTreeFormula* lFVars19 = new TTreeFormula(pExpr19,pExpr19,lTree);

  //lTree->SetBranchAddress( "jet1mprune"         , &lJP);
  //lTree->SetBranchAddress( iVar1.c_str()           , &lJT1);
  //if(iVar1 != iVar2) lTree->SetBranchAddress( iVar2.c_str()           , &lJT2); 

  int lNEvents = lTree->GetEntries();
  TFile *lOFile = new TFile("Output.root","RECREATE");
  TTree *lOTree = new TTree("DMSTree","DMSTree");
  float lMVA    = 0; lOTree->Branch("bdt_all",&lMVA ,"bdt_all/F");
  for (Long64_t i0=0; i0<lNEvents;i0++) {
    if (i0 % 10000 == 0) std::cout << "--- ... Processing event: " << double(i0)/double(lNEvents) << std::endl;
    lTree->GetEntry(i0);
    jet1QGtagComb = lFVars0->EvalInstance();
    fjet1QGtagSub1 = lFVars1->EvalInstance();
    fjet1QGtagSub2 = lFVars2->EvalInstance();
    fjet1QGtag = lFVars3->EvalInstance();
    fjet1PullAngle = lFVars4->EvalInstance();
    fjet1Pull = lFVars5->EvalInstance();
    fjet1MassTrimmed = lFVars6->EvalInstance();
    fjet1MassPruned = lFVars7->EvalInstance();
    fjet1MassSDbm1 = lFVars8->EvalInstance();
    fjet1MassSDb2 = lFVars9->EvalInstance();
    fjet1MassSDb0 = lFVars10->EvalInstance();
    fjet1QJetVol = lFVars11->EvalInstance();
    fjet1C2b2 = lFVars12->EvalInstance();
    fjet1C2b1 = lFVars13->EvalInstance();
    fjet1C2b0p5 = lFVars14->EvalInstance();
    fjet1C2b0p2 = lFVars15->EvalInstance();
    fjet1C2b0 = lFVars16->EvalInstance();
    fjet1Tau2 = lFVars17->EvalInstance();
    fjet1Tau1 = lFVars18->EvalInstance();
    tau2tau1 = lFVars19->EvalInstance();

    lMVA      = float(reader->EvaluateMVA(lJetName.c_str()));
    lOTree->Fill();
  }
  lOTree->Write();
  lOFile->Close();
  delete reader;
}

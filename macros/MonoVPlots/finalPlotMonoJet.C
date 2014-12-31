#include "TROOT.h"
#include "TSystem.h"
#include "TInterpreter.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TStyle.h"
#include "TPad.h"
#include "Math/QuantFuncMathCore.h"
#include "TMath.h"
#include "TGraphAsymmErrors.h"
#include "StandardPlotMonoJet.C"
#include <iomanip>
#include <iostream>
#include <fstream>
#include "CMS_lumi.C"

//void atributes(TH1D *histo, TString xtitle = "", TString ytitle = "Fraction");

void eraselabel(TPad *p,Double_t h){
  p->cd();
  TPad *pe = new TPad("pe","pe",0.02,0,p->GetLeftMargin(),h);	   
  pe->Draw();
  pe->SetFillColor(p->GetFillColor()); 
  pe->SetBorderMode(0);
}

void atributes(TH1D *histo, TString xtitle = "", TString ytitle = "Fraction", TString units = ""){

  histo->SetTitle("");
  //histo->SetMarkerStyle(20);
  //histo->SetMarkerSize(0.8);
  //histo->SetLineWidth(4);
  if(strcmp(units.Data(),"")==0){
    histo->GetXaxis()->SetTitle(xtitle.Data());
  } else {
    histo->GetXaxis()->SetTitle(Form("%s [%s]",xtitle.Data(),units.Data()));
  }
  histo->GetXaxis()->SetLabelFont  (   42);
  histo->GetXaxis()->SetLabelOffset(0.015);
  histo->GetXaxis()->SetLabelSize  (0.100);
  histo->GetXaxis()->SetNdivisions (  505);
  histo->GetXaxis()->SetTitleFont  (   42);
  histo->GetXaxis()->SetTitleOffset(  1.1);
  histo->GetXaxis()->SetTitleSize  (0.100);
  histo->GetXaxis()->SetTickLength (0.07 );

  histo->GetYaxis()->SetTitle(ytitle.Data());
  histo->GetYaxis()->SetLabelFont  (   42);
  histo->GetYaxis()->SetLabelOffset(0.015);
  histo->GetYaxis()->SetLabelSize  (0.120);
  histo->GetYaxis()->SetNdivisions (  505);
  histo->GetYaxis()->SetTitleFont  (   42);
  histo->GetYaxis()->SetTitleOffset(  0.5);
  histo->GetYaxis()->SetTitleSize  (0.120);
  histo->GetYaxis()->SetTickLength (0.03 );

  histo->SetLineColor  (kBlack);
  histo->SetMarkerStyle(kFullCircle);
}

void finalPlotMonoJet(int nsel = 0, int ReBin = 1, TString XTitle = "N_{jets}", TString units = "", 
                      TString plotName = "BDT_loose_metRaw.root", TString outputName = "njets",
                      bool isLogY = false, double lumi = 19.7, bool isBlinded = true, bool scaleToData = false,
                      bool isVarBins = false) {

  gInterpreter->ExecuteMacro("GoodStyle.C");
  gStyle->SetOptStat(0);
  TFile* file = new TFile(plotName.Data(), "read");
  if(!file->IsOpen()) return;

  StandardPlot myPlot;
  myPlot.setBlinded(isBlinded);
  myPlot.setLabel(XTitle.Data());
  myPlot.setVarBins(isVarBins);
  if     (lumi ==    4.9) myPlot.setLumiLabel("4.9 fb^{-1} (7 TeV)");
  else if(lumi ==   19.7) myPlot.setLumiLabel("19.7 fb^{-1} (8 TeV)");
  else if(lumi ==   24.4) myPlot.setLumiLabel("4.9 fb^{-1} (7 TeV) + 19.7 fb^{-1} (8 TeV)");
  else                    myPlot.setLumiLabel(""); 
  myPlot.setUnits(units.Data());

  TH1D* hTop     = (TH1D*)file->Get("Top");
  TH1D* hWjets   = (TH1D*)file->Get("W+jets");
  TH1D* hZjets   = (TH1D*)file->Get("Z+jets");
  TH1D* hDiboson = (TH1D*)file->Get("Diboson");
  TH1D* hOther   = (TH1D*)file->Get("Others");
  TH1D* hData    = (TH1D*)file->Get("htempdata_0");
  TH1D* hGjets = 0;
  if (nsel == 3)
    hGjets = (TH1D*)file->Get("G+jets"); 
  
  TH1D* hTotal = (TH1D*) hTop->Clone();
  hTotal->Reset();
  TH1D* hTotalDivision = (TH1D*) hTop->Clone();
  hTotalDivision->Reset();
  TH1D* hDataDivision = (TH1D*) hTop->Clone();
  hDataDivision->Reset();
  TH1D* hRatio = (TH1D*) hTop->Clone();
  hRatio->Reset();
  
  TH1D* hHiggs = 0;
  if (nsel == 0) {
    hHiggs = (TH1D*)file->Get("WH_125");
    hHiggs->Add((TH1D*)file->Get("ZH_125"));
    hHiggs->Add((TH1D*)file->Get("ggH_125"));
  }

  int binToRemove = -1; //hTop->GetNbinsX();
  hTop     ->SetBinContent(binToRemove,0);
  hWjets   ->SetBinContent(binToRemove,0);
  hZjets   ->SetBinContent(binToRemove,0);
  hDiboson ->SetBinContent(binToRemove,0);
  hOther   ->SetBinContent(binToRemove,0);
  hData    ->SetBinContent(binToRemove,0);
  if (nsel == 3) 
    hGjets ->SetBinContent(binToRemove,0);


  Double_t totalBkg = hTop->GetSumOfWeights()+hWjets->GetSumOfWeights()+hZjets->GetSumOfWeights()+
                      hDiboson->GetSumOfWeights()+hOther->GetSumOfWeights();
  if (nsel == 3) 
    totalBkg += hGjets->GetSumOfWeights();

  double scale = hData->GetSumOfWeights()/totalBkg;
  if (scaleToData) {
    printf("(normalizing MC to data: %.2f) - %.2f\n",scale,totalBkg);
    hTop    ->Scale(scale);
    hWjets  ->Scale(scale);
    hZjets  ->Scale(scale);
    hDiboson->Scale(scale);
    hOther  ->Scale(scale);
    if (nsel == 3)
      hGjets ->Scale(scale);
  }

  // Construct the total background
  hTotal->Add(hOther  );
  hTotal->Add(hDiboson);
  hTotal->Add(hTop    );
  hTotal->Add(hWjets  );
  hTotal->Add(hZjets  );
  if (nsel == 3)
    hTotal->Add(hGjets);

  Int_t theNselPlot = nsel;
  
  if(nsel == 0 && hHiggs->GetSumOfWeights() > 0){
    if (scaleToData) 
      hHiggs->Scale(totalBkg/hHiggs->GetSumOfWeights()*scale);
    myPlot.setMCHist(iHiggs,(TH1D*)hHiggs->Clone("hHiggs"));
  }
  if     (nsel == 0){ // default
    myPlot.setMCHist(iOther,  (TH1D*)hOther  ->Clone("hOther"));
    myPlot.setMCHist(iDiboson,(TH1D*)hDiboson->Clone("hDiboson"));
    myPlot.setMCHist(iTop,    (TH1D*)hTop    ->Clone("hTop"));
    myPlot.setMCHist(iWjets,  (TH1D*)hWjets  ->Clone("hWjets"));
    myPlot.setMCHist(iZjets,  (TH1D*)hZjets  ->Clone("hZjets")); 
    hData->Reset();
    myPlot.setDataHist((TH1D*)hData->Clone("data"));
  }
  else if (nsel == 1){ // zll
    myPlot.setMCHist(iDiboson,(TH1D*)hDiboson->Clone("hDiboson"));
    myPlot.setMCHist(iTop,    (TH1D*)hTop    ->Clone("hTop"));
    myPlot.setMCHist(iWjets,  (TH1D*)hWjets  ->Clone("hWjets"));
    myPlot.setMCHist(iZjets,  (TH1D*)hZjets  ->Clone("hZjets")); 
    myPlot.setDataHist((TH1D*)hData->Clone("data"));
  }
  else if (nsel == 2){ // wlv
    myPlot.setMCHist(iOther,  (TH1D*)hOther  ->Clone("hOther"));
    myPlot.setMCHist(iDiboson,(TH1D*)hDiboson->Clone("hDiboson"));
    myPlot.setMCHist(iTop,    (TH1D*)hTop    ->Clone("hTop"));
    myPlot.setMCHist(iWjets,  (TH1D*)hWjets  ->Clone("hWjets"));
    myPlot.setMCHist(iZjets,  (TH1D*)hZjets  ->Clone("hZjets")); 
    myPlot.setDataHist((TH1D*)hData->Clone("data"));
  }
  else if (nsel == 3){ // pj
    hOther->Add(hDiboson);
    hOther->Add(hTop    );
    hOther->Add(hWjets  );
    hOther->Add(hZjets  );
    myPlot.setMCHist(iGjets,  (TH1D*)hGjets  ->Clone("hGjets"));
    myPlot.setMCHist(iOther,  (TH1D*)hOther  ->Clone("hOther"));
    myPlot.setDataHist((TH1D*)hData->Clone("data"));
  }
  else assert(0);
  
  myPlot.setNsel(theNselPlot);

  TCanvas* c1 = new TCanvas("c1", "c1",5,5,500,500);
  c1->SetBottomMargin(0.1);
  c1->cd();
  TPad *pad1 = new TPad("pad1", "pad1",0.00,0.30,1.00,1.00);
  TPad *pad2 = new TPad("pad2", "pad2",0.00,0.00,1.00,0.30);
  pad1->SetBottomMargin(0);
  pad2->SetTopMargin(0);
  pad2->SetBottomMargin(0.3);
  pad1->Draw();
  pad2->Draw();

  pad1->cd();
  if(isLogY == true) c1->SetLogy();
  if(isLogY == true) pad1->SetLogy();
  myPlot.Draw(ReBin);  // Can pass a rebin 
  CMS_lumi( pad1, 2, 12 );

  pad2->cd();
  hDataDivision ->Add(hData );
  hTotalDivision->Add(hTotal);
  hDataDivision ->Rebin(ReBin);
  hTotalDivision->Rebin(ReBin);
  hRatio        ->Rebin(ReBin);
  // Compute pull with systematic uncertainties
  double systBck = 0.06;
  double poissErrNoEntries = ROOT::Math::gamma_quantile_c(1-0.6827,1,1.);
  for(int i=1; i<=hData->GetNbinsX(); i++){
    double errorSq = 0;
    // Add minimal data entries in case of no observed events
    if (!(hDataDivision->GetBinError(i) > 0))
      errorSq += poissErrNoEntries*poissErrNoEntries;
    double pull = 0.0;
    errorSq += hDataDivision->GetBinError(i)*hDataDivision->GetBinError(i)+hTotalDivision->GetBinError(i)*hTotalDivision->GetBinError(i);
    errorSq += hTotalDivision->GetBinContent(i)*hTotalDivision->GetBinContent(i)*systBck*systBck;
    pull = (hDataDivision->GetBinContent(i)-hTotalDivision->GetBinContent(i))/sqrt(errorSq);    
    hRatio->SetBinContent(i,pull);
    hRatio->SetBinError(i,1.0);
    if (isBlinded) {
      hRatio->SetBinContent(i,0.0);
      hRatio->SetBinError(i,0.0);
    }    
  }
  atributes(hRatio,XTitle.Data(),"Pull",units.Data());
  hRatio->Draw("e");

  // Draw a line throgh y=0
  TLine* baseline = new TLine(hRatio->GetXaxis()->GetXmin(), 0.,
                              hRatio->GetXaxis()->GetXmax(), 0.);
  baseline->SetLineStyle(kDashed);
  baseline->Draw();
  // Set the y-axis range symmetric around y=0
  Double_t dy = TMath::Max(TMath::Abs(hRatio->GetMaximum()),
                           TMath::Abs(hRatio->GetMinimum())) + 1.5;
  hRatio->GetYaxis()->SetRangeUser(-dy, +dy);
  hRatio->GetYaxis()->CenterTitle();
  eraselabel(pad1,hData->GetXaxis()->GetLabelSize());

  char CommandToExec[300];
  sprintf(CommandToExec,"mkdir -p plots");
  gSystem->Exec(CommandToExec);  

  if(strcmp(outputName.Data(),"")!=0){
    char myOutputFile[300];
    sprintf(myOutputFile,"plots/%s.eps",outputName.Data());
    c1->SaveAs(myOutputFile);
    sprintf(myOutputFile,"plots/%s.png",outputName.Data());
    c1->SaveAs(myOutputFile);
    sprintf(myOutputFile,"plots/%s.pdf",outputName.Data());
    c1->SaveAs(myOutputFile);
    sprintf(myOutputFile,"plots/%s.root",outputName.Data());
    c1->SaveAs(myOutputFile);
  }

}

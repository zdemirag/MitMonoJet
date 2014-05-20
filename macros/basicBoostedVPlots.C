//--------------------------------------------------------------------------------------------------
// Basic plots for the boostedV analysis.
//
// Authors: C.Paus                                                                        (Mar 2014)
//--------------------------------------------------------------------------------------------------
#include <TSystem.h>
#include "MitPlots/Style/interface/MitStyle.h"
#include "MitPlots/Plot/interface/PlotTask.h"
using namespace std;
using namespace mithep;
//==================================================================================================
void basicBoostedVPlots(double lumi = 19600.0)
{
  // Create the basic plots of the boosted V analysis to check the key N subjettiness variables


  // set the folder containing the input ntuples properly
  // here you can change the plot sources, these are the defaults
  // gSystem->Setenv("MIT_PROD_CFG","boostedv");
  // gSystem->Setenv("MIT_ANA_HIST","/home/$USER/cms/hist/boostedv/merged");

  // setup graphics stuff before starting
  MitStyle::Init();

  // plot from TTree
  TString nTuple   = "DMSTree";
  printf("\n Looking at tree (%s)\n",nTuple.Data());

  // predefined cuts
  TString orC(" || ");
  TString andC(" && ");
  TString leptonCuts("(lep1.Pt() > 30 && abs(lep1.Eta()) < 2.1 && abs(lid1) > 12)");
  TString jetCuts("(tjet.Pt() > 300 && abs(tjet.Eta()) < 2.4 && nbjets > 1)");
  TString nSubCuts("(tjetTau2/tjetTau1 < 0.5)");
  TString ECFCuts("(tjetC2b0 > 100.)");
  TString eventWeight("(puweight)*(");

  // our plot task
  TString variable, cuts;
  PlotTask *plotTask = 0;

  // set plot config properly
  gSystem->Setenv("MIT_ANA_CFG","boostedv-plots-top");

  // raw met with basic cuts
  plotTask = new PlotTask(0,lumi);
  variable = "metRaw";
  cuts     = eventWeight+jetCuts+andC+leptonCuts+")";
  plotTask->SetHistRanges(0.0,500.,0.,0.);
  plotTask->SetNBins(25);
  plotTask->SetAxisTitles("raw PFMET [GeV]","Number of Events");
  plotTask->SetPngFileName("rawMet.png");
  plotTask->Plot(Stacked,nTuple,variable,cuts);
  delete plotTask;

  // met mva with basic cuts
  plotTask = new PlotTask(0,lumi);
  variable = "mvamet";
  cuts     = eventWeight+jetCuts+andC+leptonCuts+")";
  plotTask->SetHistRanges(0.0,500.,0.,0.);
  plotTask->SetNBins(25);
  plotTask->SetAxisTitles("MVA MET [GeV]","Number of Events");
  plotTask->SetPngFileName("mvaMet.png");
  plotTask->Plot(Stacked,nTuple,variable,cuts);
  delete plotTask;

  // raw met with basic cuts
  plotTask = new PlotTask(0,lumi);
  variable = "metRaw";
  cuts     = eventWeight+jetCuts+andC+leptonCuts+")";
  plotTask->SetHistRanges(0.0,500.,0.,0.);
  plotTask->SetLogy(kTRUE);
  plotTask->SetNBins(25);
  plotTask->SetAxisTitles("raw PFMET [GeV]","Number of Events");
  plotTask->SetPngFileName("rawMet_log.png");
  plotTask->Plot(Stacked,nTuple,variable,cuts);
  delete plotTask;

  // met mva with basic cuts
  plotTask = new PlotTask(0,lumi);
  variable = "mvamet";
  cuts     = eventWeight+jetCuts+andC+leptonCuts+")";
  plotTask->SetHistRanges(0.0,500.,0.,0.);
  plotTask->SetLogy(kTRUE);
  plotTask->SetNBins(25);
  plotTask->SetAxisTitles("MVA MET [GeV]","Number of Events");
  plotTask->SetPngFileName("mvaMet_log.png");
  plotTask->Plot(Stacked,nTuple,variable,cuts);
  delete plotTask;

  // lep pt with basic cuts
  plotTask = new PlotTask(0,lumi);
  variable = "lep1.Pt()";
  cuts     = eventWeight+jetCuts+andC+leptonCuts+")";
  plotTask->SetHistRanges(30,300,0.,0.);
  plotTask->SetNBins(27);
  plotTask->SetAxisTitles("lepton p_{T} [GeV]","Number of Events");
  plotTask->SetPngFileName("lepPt.png");
  plotTask->Plot(Stacked,nTuple,variable,cuts);
  delete plotTask;

  // lep eta with basic cuts
  plotTask = new PlotTask(0,lumi);
  variable = "lep1.Eta()";
  cuts     = eventWeight+jetCuts+andC+leptonCuts+")";
  plotTask->SetHistRanges(-2.1,2.1,0.,0.);
  plotTask->SetNBins(30);
  plotTask->SetAxisTitles("lepton #eta","Number of Events");
  plotTask->SetPngFileName("lepEta.png");
  plotTask->Plot(Stacked,nTuple,variable,cuts);
  delete plotTask;

  // lep-MET dphi with basic cuts
  plotTask = new PlotTask(0,lumi);
  variable = "TMath::ACos(cos(lep1.Phi()-metRawPhi))";
  cuts     = eventWeight+jetCuts+andC+leptonCuts+")";
  plotTask->SetHistRanges(0,3.2,0.,0.);
  plotTask->SetNBins(30);
  plotTask->SetAxisTitles("#Delta#phi (lep-PFMET)","Number of Events");
  plotTask->SetPngFileName("lepRawMetDphi.png");
  plotTask->Plot(Stacked,nTuple,variable,cuts);
  delete plotTask;

  // lep-MET dphi with basic cuts
  plotTask = new PlotTask(0,lumi);
  variable = "TMath::ACos(cos(lep1.Phi()-mvametPhi))";
  cuts     = eventWeight+jetCuts+andC+leptonCuts+")";
  plotTask->SetHistRanges(0,3.2,0.,0.);
  plotTask->SetNBins(30);
  plotTask->SetAxisTitles("#Delta#phi (lep-MVAMET)","Number of Events");
  plotTask->SetPngFileName("lepMvaMetDphi.png");
  plotTask->Plot(Stacked,nTuple,variable,cuts);
  delete plotTask;

  // jet pt with basic cuts
  plotTask = new PlotTask(0,lumi);
  variable = "tjet.Pt()";
  cuts     = eventWeight+jetCuts+andC+leptonCuts+")";
  plotTask->SetHistRanges(300.0,1000.,0.,0.);
  plotTask->SetNBins(50);
  plotTask->SetAxisTitles("jet p_{T} [GeV]","Number of Events");
  plotTask->SetPngFileName("tjPt.png");
  plotTask->Plot(Stacked,nTuple,variable,cuts);
  delete plotTask;

  // jet b-tag with basic cuts
  plotTask = new PlotTask(0,lumi);
  variable = "tjet.Eta()";
  cuts     = eventWeight+jetCuts+andC+leptonCuts+")";
  plotTask->SetHistRanges(-2.4,2.4,0.,0.);
  plotTask->SetNBins(30);
  plotTask->SetAxisTitles("jet #eta","Number of Events");
  plotTask->SetPngFileName("tjEta.png");
  plotTask->Plot(Stacked,nTuple,variable,cuts);
  delete plotTask;

  // jet eta with basic cuts
  plotTask = new PlotTask(0,lumi);
  variable = "tjetBtag";
  cuts     = eventWeight+jetCuts+andC+leptonCuts+")";
  plotTask->SetHistRanges(0,1,0.,0.);
  plotTask->SetNBins(50);
  plotTask->SetAxisTitles("jet b-tag (CSV)","Number of Events");
  plotTask->SetPngFileName("tjBtag.png");
  plotTask->Plot(Stacked,nTuple,variable,cuts);
  delete plotTask;

  // Jet multiplicity with basic cuts
  plotTask = new PlotTask(0,lumi);
  variable = "njets";
  cuts     = eventWeight+jetCuts+andC+leptonCuts+")";
  plotTask->SetHistRanges(0.0,10.,0.,0.);
  plotTask->SetNBins(10);
  plotTask->SetAxisTitles("jet multiplicity","Number of Events");
  plotTask->SetPngFileName("njets.png");
  plotTask->Plot(Stacked,nTuple,variable,cuts);
  delete plotTask;

  // b-Jets multiplicity with basic cuts
  plotTask = new PlotTask(0,lumi);
  variable = "nbjets";
  cuts     = eventWeight+jetCuts+andC+leptonCuts+")";
  plotTask->SetHistRanges(0.0,5.,0.,0.);
  plotTask->SetNBins(5);
  plotTask->SetAxisTitles("b-jet multiplicity","Number of Events");
  plotTask->SetPngFileName("nbjets.png");
  plotTask->Plot(Stacked,nTuple,variable,cuts);
  delete plotTask;

  // Mass with basic cuts
  plotTask = new PlotTask(0,lumi);
  variable = "tjet.M()";
  cuts     = eventWeight+jetCuts+andC+leptonCuts+")";
  plotTask->SetHistRanges(0.0,200.,0.,0.);
  plotTask->SetNBins(50);
  plotTask->SetAxisTitles("jet mass [GeV]","Number of Events");
  plotTask->SetPngFileName("tjMass.png");
  plotTask->Plot(Stacked,nTuple,variable,cuts);
  delete plotTask;

  // Groomed Mass with basic cuts
  plotTask = new PlotTask(0,lumi);
  variable = "tjetGroomed.M()";
  cuts     = eventWeight+jetCuts+andC+leptonCuts+")";
  plotTask->SetHistRanges(0.0,200.,0.,0.);
  plotTask->SetNBins(50);
  plotTask->SetAxisTitles("groomed jet mass [GeV]","Number of Events");
  plotTask->SetPngFileName("tjGroomedMass.png");
  plotTask->Plot(Stacked,nTuple,variable,cuts);
  delete plotTask;

  // Subjettiness with basic cuts
  plotTask = new PlotTask(0,lumi);
  variable = "tjetTau2/tjetTau1";
  cuts     = eventWeight+jetCuts+andC+leptonCuts+")";
  plotTask->SetHistRanges(0.0,1.,0.,0.);
  plotTask->SetNBins(50);
  plotTask->SetAxisTitles("jet #tau2/#tau1","Number of Events");
  plotTask->SetPngFileName("tjT2overT1.png");
  plotTask->Plot(Stacked,nTuple,variable,cuts);
  delete plotTask;

  // ECF with basic cuts
  plotTask = new PlotTask(0,lumi);
  variable = "tjetC2b0";
  cuts     = eventWeight+jetCuts+andC+leptonCuts+")";
  plotTask->SetHistRanges(0.0,200.,0.,0.);
  plotTask->SetNBins(50);
  plotTask->SetAxisTitles("jet C_{2} (#beta = 0)","Number of Events");
  plotTask->SetPngFileName("tjC2b0.png");
  plotTask->Plot(Stacked,nTuple,variable,cuts);
  delete plotTask;

  // Mass after subjettiness cut
  plotTask = new PlotTask(0,lumi);
  variable = "tjet.M()";
  cuts     = eventWeight+jetCuts+andC+leptonCuts+andC+nSubCuts+")";
  plotTask->SetHistRanges(0.0,200.,0.,0.);
  plotTask->SetNBins(50);
  plotTask->SetAxisTitles("jet mass [GeV]","Number of Events");
  plotTask->SetPngFileName("tjMass-nsubCut.png");
  plotTask->Plot(Stacked,nTuple,variable,cuts);
  delete plotTask;

  // Groomed Mass after subjettiness cut
  plotTask = new PlotTask(0,lumi);
  variable = "tjetGroomed.M()";
  cuts     = eventWeight+jetCuts+andC+leptonCuts+andC+nSubCuts+")";
  plotTask->SetHistRanges(0.0,200.,0.,0.);
  plotTask->SetNBins(50);
  plotTask->SetAxisTitles("groomed jet mass [GeV]","Number of Events");
  plotTask->SetPngFileName("tjGroomedMass-nsubCut.png");
  plotTask->Plot(Stacked,nTuple,variable,cuts);
  delete plotTask;

  // Mass after ECF cut
  plotTask = new PlotTask(0,lumi);
  variable = "tjet.M()";
  cuts     = eventWeight+jetCuts+andC+leptonCuts+andC+ECFCuts+")";
  plotTask->SetHistRanges(0.0,200.,0.,0.);
  plotTask->SetNBins(50);
  plotTask->SetAxisTitles("jet mass [GeV]","Number of Events");
  plotTask->SetPngFileName("tjMass-ECFCut.png");
  plotTask->Plot(Stacked,nTuple,variable,cuts);
  delete plotTask;

  // Groomed Mass after ECF cut
  plotTask = new PlotTask(0,lumi);
  variable = "tjetGroomed.M()";
  cuts     = eventWeight+jetCuts+andC+leptonCuts+andC+ECFCuts+")";
  plotTask->SetHistRanges(0.0,200.,0.,0.);
  plotTask->SetNBins(50);
  plotTask->SetAxisTitles("groomed jet mass [GeV]","Number of Events");
  plotTask->SetPngFileName("tjGroomedMass-ECFCut.png");
  plotTask->Plot(Stacked,nTuple,variable,cuts);
  delete plotTask;

  return;
}

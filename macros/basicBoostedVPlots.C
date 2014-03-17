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
void basicBoostedVPlots(double lumi = 4500.0)
{
  // setup graphics stuff before starting
  MitStyle::Init();

  // plot from TTree
  TString nTuple   = "MJetTree";
  printf("\n Looking at tree (%s)\n",nTuple.Data());

  // predefined cuts
  TString orC(" || ");
  TString andC(" && ");
  TString leptonCuts("nLep+nTaus<1");
  TString leptonMcCuts("genLep1Pid==0");
  TString basicCuts =
    TString("genJet1Pt>300 && genJet1M>65. && genJet1M<95. &&")+
    TString("genJet1Eta<1.3 && genJet1Eta>-1.3");

  // our plot task
  TString variable, cuts;
  PlotTask *plotTask = 0;


  // set plot config properly
  gSystem->Setenv("MIT_ANA_CFG","boostedv-plots-w");

  // Tau2/Tau1 - basic cuts
  plotTask = new PlotTask(0,lumi);
  variable = "genJet1Tau2/genJet1Tau1";
  cuts     = basicCuts;
  plotTask->SetHistRanges(0.0,1.0,0.,0.);
  plotTask->SetAxisTitles("#tau_{2}/#tau_{1}","Number of Events");
  plotTask->SetPngFileName("t2OverT1.png");
  plotTask->Plot(Normalized,nTuple,variable,cuts);
  delete plotTask;
  
  // Tau2/Tau1 - basic cuts and reco leptons rejected
  plotTask = new PlotTask(0,lumi);
  variable = "genJet1Tau2/genJet1Tau1";
  cuts     = basicCuts+andC+leptonCuts;
  plotTask->SetHistRanges(0.0,1.0,0.,0.);
  plotTask->SetAxisTitles("#tau_{2}/#tau_{1}","Number of Events");
  plotTask->SetPngFileName("t2OverT1Lep.png");
  plotTask->Plot(Normalized,nTuple,variable,cuts);
  delete plotTask;

  // Tau2/Tau1 - basic cuts and MC leptons rejected
  plotTask = new PlotTask(0,lumi);
  variable = "genJet1Tau2/genJet1Tau1";
  cuts     = basicCuts+andC+leptonMcCuts;
  plotTask->SetHistRanges(0.0,1.0,0.,1.);
  plotTask->SetAxisTitles("#tau_{2}/#tau_{1}","Number of Events");
  plotTask->SetPngFileName("t2OverT1LepMc.png");
  plotTask->Plot(Normalized,nTuple,variable,cuts);
  delete plotTask;

  // set plot config properly
  gSystem->Setenv("MIT_ANA_CFG","boostedv-plots-top");

  basicCuts =
    TString("genJet1Pt>300 && genJet1M>145. && genJet1M<205. &&")+
    TString("genJet1Eta<1.3 && genJet1Eta>-1.3");

  // Tau3/Tau2 - basic cuts
  plotTask = new PlotTask(0,lumi);
  variable = "genJet1Tau3/genJet1Tau2";
  cuts     = basicCuts;
  plotTask->SetHistRanges(0.0,1.0,0.,0.);
  plotTask->SetAxisTitles("#tau_{3}/#tau_{2}","Number of Events");
  plotTask->SetPngFileName("t3OverT2.png");
  plotTask->Plot(Normalized,nTuple,variable,cuts);
  delete plotTask;
  
  // Tau2/Tau1 - basic cuts and reco leptons rejected
  plotTask = new PlotTask(0,lumi);
  variable = "genJet1Tau3/genJet1Tau2";
  cuts     = basicCuts+andC+leptonCuts;
  plotTask->SetHistRanges(0.0,1.0,0.,0.);
  plotTask->SetAxisTitles("#tau_{3}/#tau_{2}","Number of Events");
  plotTask->SetPngFileName("t3OverT2Lep.png");
  plotTask->Plot(Normalized,nTuple,variable,cuts);
  delete plotTask;

  // Tau2/Tau1 - basic cuts and MC leptons rejected
  plotTask = new PlotTask(0,lumi);
  variable = "genJet1Tau3/genJet1Tau2";
  cuts     = basicCuts+andC+leptonMcCuts;
  plotTask->SetHistRanges(0.0,1.0,0.,1.);
  plotTask->SetAxisTitles("#tau_{3}/#tau_{2}","Number of Events");
  plotTask->SetPngFileName("t3OverT2LepMc.png");
  plotTask->Plot(Normalized,nTuple,variable,cuts);
  delete plotTask;

  return;
}

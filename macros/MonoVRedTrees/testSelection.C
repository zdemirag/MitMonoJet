//--------------------------------------------------------------------------------------------------
// Basic plots for the boostedV analysis.
//
// Original author: C.Paus                                                                (Mar 2014)
// Additional authors: L.Di Matteo, B.Allen                                               (Aug 2014)
//--------------------------------------------------------------------------------------------------
#include <TSystem.h>
#include "MitCommon/MathTools/interface/MathUtils.h"
#include "MitPlots/Style/interface/MitStyle.h"
#include "MitPlots/Plot/interface/PlotTask.h"
using namespace std;
using namespace mithep;

//==================================================================================================
void testSelection(double lumi = 19700.0)
{
  gSystem->Setenv("MIT_ANA_CFG","boostedv-ana-test");
  gSystem->Setenv("MIT_ANA_HIST","/mnt/hscratch/dimatteo/boostedv-v9/merged/");
 
  // plot from TTree
  TString nTuple   = "DMSTree";
  printf("\n Looking at tree (%s)\n",nTuple.Data());

  //TString cut("(puweight)*(((trigger & (1<<0)) || (trigger & (1<<1))) && ((jet1.Pt() > 110 && abs(jet1.eta()) < 2.5 && jet1CHF > 0.2 && jet1NHF < 0.7 && jet1NEMF < 0.7)) && (njets == 1 || (njets == 2 && (TMath::ACos(TMath::Cos(jet2.phi() - jet1.phi())) < 2.0))) && (ntaus == 0) && (nphotons == 0) && (metFiltersWord == 511 || metFiltersWord == 1023) && (fjet1.Pt() > 250 && abs(fjet1.Eta()) < 2.5) && (fjet1Tau2/fjet1Tau1 < 0.5 && 60 < fjet1MassPruned && fjet1MassPruned < 110) && (metRaw > 250 && metRaw < 1000) && (preselWord & (1<<3)) && (nlep == 0) && (metRaw > 200 && metRaw < 1000))");
  //TString cut("(puweight)*(((trigger & (1<<0)) || (trigger & (1<<1))) && ((jet1.Pt() > 110 && abs(jet1.eta()) < 2.5 && jet1CHF > 0.2 && jet1NHF < 0.7 && jet1NEMF < 0.7)) && (njets == 1 || (njets == 2 && (TMath::ACos(TMath::Cos(jet2.phi() - jet1.phi())) < 2.0))) && (ntaus == 0) && (nphotons == 0) && (metFiltersWord == 511 || metFiltersWord == 1023) && (fjet1.Pt() > 250 && abs(fjet1.Eta()) < 2.5) && (fjet1Tau2/fjet1Tau1 < 0.5 && 60 < fjet1MassPruned && fjet1MassPruned < 110) && (preselWord & (1<<2)) && (nlep == 2 && (abs(lid1)>=130 && abs(lid2)>=13) && abs(lep1.Eta()) < 2.1 && abs(lep2.Eta()) < 2.1) && (60 < mll && mll < 120) && (met > 250 && met < 1000))");
  //TString cut("(puweight)*(((trigger & (1<<0)) || (trigger & (1<<1))) && ((jet1.Pt() > 110 && abs(jet1.eta()) < 2.5 && jet1CHF > 0.2 && jet1NHF < 0.7 && jet1NEMF < 0.7)) && (njets == 1 || (njets == 2 && (TMath::ACos(TMath::Cos(jet2.phi() - jet1.phi())) < 2.0))) && (ntaus == 0) && (nphotons == 0) && (metFiltersWord == 511 || metFiltersWord == 1023) && (fjet1.Pt() > 250 && abs(fjet1.Eta()) < 2.5) && (fjet1Tau2/fjet1Tau1 < 0.5 && 60 < fjet1MassPruned && fjet1MassPruned < 110) && (preselWord & (1<<1)) && (nlep == 1 && abs(lid1) >= 1300 && lep1.Pt() > 20 && abs(lep1.Eta()) < 2.1) && (40 < mt && mt < 200) && (met > 250 && met < 1000))");
  //TString cut("(puweight)*((trigger & (1<<3)) && (preselWord & (1<<5)) && (nlep == 0) && ((jet1.Pt() > 110 && abs(jet1.eta()) < 2.5 && jet1CHF > 0.2 && jet1NHF < 0.7 && jet1NEMF < 0.7)) && (njets == 1 || (njets == 2 && (TMath::ACos(TMath::Cos(jet2.phi() - jet1.phi())) < 2.0))) && (ntaus == 0) && (metFiltersWord == 511 || metFiltersWord == 1023) && (pho1.Pt() > 160 && abs(pho1.Eta()) < 2.5) && (met > 250 && met < 1000) && (fjet1.Pt() > 250 && abs(fjet1.Eta()) < 2.5) && (fjet1Tau2/fjet1Tau1 < 0.5 && 60 < fjet1MassPruned && fjet1MassPruned < 110))");
  
  //TString cut("(puweight)*(((trigger & (1<<0)) || (trigger & (1<<1))) && ((jet1.Pt() > 110 && abs(jet1.eta()) < 2.5 && jet1CHF > 0.2 && jet1NHF < 0.7 && jet1NEMF < 0.7)) && (njets == 1 || (njets == 2 && (TMath::ACos(TMath::Cos(jet2.phi() - jet1.phi())) < 2.0))) && (ntaus == 0) && (nphotons == 0) && (metFiltersWord == 511 || metFiltersWord == 1023) && (jet1.Pt() > 150 && abs(jet1.Eta()) < 2.0) && !((fjet1.Pt() > 250 && abs(fjet1.Eta()) < 2.5) && (fjet1Tau2/fjet1Tau1 < 0.5 && 60 < fjet1MassPruned && fjet1MassPruned < 110) && (met > 250)) && (preselWord & (1<<3)) && (nlep == 0) && (metRaw > 200 && metRaw < 1000))");
  //TString cut("(puweight)*(((trigger & (1<<0)) || (trigger & (1<<1))) && ((jet1.Pt() > 110 && abs(jet1.eta()) < 2.5 && jet1CHF > 0.2 && jet1NHF < 0.7 && jet1NEMF < 0.7)) && (njets == 1 || (njets == 2 && (TMath::ACos(TMath::Cos(jet2.phi() - jet1.phi())) < 2.0))) && (ntaus == 0) && (nphotons == 0) && (metFiltersWord == 511 || metFiltersWord == 1023) && (jet1.Pt() > 150 && abs(jet1.Eta()) < 2.0) && !((fjet1.Pt() > 250 && abs(fjet1.Eta()) < 2.5) && (fjet1Tau2/fjet1Tau1 < 0.5 && 60 < fjet1MassPruned && fjet1MassPruned < 110) && (met > 250)) && (preselWord & (1<<2)) && (nlep == 2 && (abs(lid1)>=130 && abs(lid2)>=13) && abs(lep1.Eta()) < 2.1 && abs(lep2.Eta()) < 2.1) && (60 < mll && mll < 120) && (met > 200 && met < 1000))");
  //TString cut("(puweight)*(((trigger & (1<<0)) || (trigger & (1<<1))) && ((jet1.Pt() > 110 && abs(jet1.eta()) < 2.5 && jet1CHF > 0.2 && jet1NHF < 0.7 && jet1NEMF < 0.7)) && (njets == 1 || (njets == 2 && (TMath::ACos(TMath::Cos(jet2.phi() - jet1.phi())) < 2.0))) && (ntaus == 0) && (nphotons == 0) && (metFiltersWord == 511 || metFiltersWord == 1023) && (jet1.Pt() > 150 && abs(jet1.Eta()) < 2.0) && !((fjet1.Pt() > 250 && abs(fjet1.Eta()) < 2.5) && (fjet1Tau2/fjet1Tau1 < 0.5 && 60 < fjet1MassPruned && fjet1MassPruned < 110) && (met > 250)) && (preselWord & (1<<1)) && (nlep == 1 && abs(lid1) >= 1300 && lep1.Pt() > 20 && abs(lep1.Eta()) < 2.1) && (40 < mt && mt < 200) && (met > 200 && met < 1000))");
  TString cut("(puweight)*((trigger & (1<<3)) && (preselWord & (1<<5)) && (nlep == 0) && ((jet1.Pt() > 110 && abs(jet1.eta()) < 2.5 && jet1CHF > 0.2 && jet1NHF < 0.7 && jet1NEMF < 0.7)) && (njets == 1 || (njets == 2 && (TMath::ACos(TMath::Cos(jet2.phi() - jet1.phi())) < 2.0))) && (ntaus == 0) && (metFiltersWord == 511 || metFiltersWord == 1023) && (pho1.Pt() > 160 && abs(pho1.Eta()) < 2.5) && (met > 200 && met < 1000) && (jet1.Pt() > 150 && abs(jet1.Eta()) < 2.0) && !((fjet1.Pt() > 250 && abs(fjet1.Eta()) < 2.5) && (fjet1Tau2/fjet1Tau1 < 0.5 && 60 < fjet1MassPruned && fjet1MassPruned < 110) && (met > 250)))");
  
  // plot metRaw
  PlotTask *plotTask = 0;
  plotTask = new PlotTask(0,lumi);
  plotTask->SetVarBins(false);
  plotTask->SetNormToWidth(false);
  plotTask->SetHistRanges(0.,2.,0.,0.);
  plotTask->SetNBins(2);
  plotTask->SetAxisTitles("isData","Number of Events");
  plotTask->SetPngFileName("/tmp/dummy.png");
  plotTask->Plot(Stacked,nTuple,"isData",cut,"test");
  delete plotTask;
  
  return;
}

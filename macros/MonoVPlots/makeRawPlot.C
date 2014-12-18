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
void makeRawPlot(double lumi = 19700.0, int mode = 0, TString variable = "metRaw", int icut = 0, int nBins = 60, float min = 200, float max = 800, bool varbins = false)
{
  // Create raw plots for further analysis (limits) or plotting

  // set the folder containing the input ntuples properly
  // here you can change the plot sources, i.e. the list of samples to be processed
  // and the location of the input flat trees
  gSystem->Setenv("MIT_ANA_CFG","boostedv-plots");
  gSystem->Setenv("MIT_ANA_HIST","/scratch4/dimatteo/cms/hist/boostedv-v9/merged-test/");
  //gSystem->Setenv("MIT_ANA_HIST","/scratch/dabercro/boostedv-v9/WithBDT/");

  // setup graphics stuff before starting
  MitStyle::Init();

  // plot from TTree
  TString nTuple   = "DMSTree";
  printf("\n Looking at tree (%s)\n",nTuple.Data());

  // predefined cuts
  TString orC(" || ");
  TString andC(" && ");
    
  // trigger and preselection cuts
  TString eventWeight = "(puweight)*("; //remember to put the parenthesis back in-
  TString monoJetTrigger = "((trigger & (1<<0)) || (trigger & (1<<1)))";
  TString singleMuTrigger = "(trigger & (1<<2))";
  TString singlePhTrigger = "(trigger & (1<<3))";
  TString metFiltersCut = "(metFiltersWord == 511 || metFiltersWord == 1023)"; //all MET filters applied but the loose beam halo
  TString preselCutsZnunu = "(preselWord & (1<<3))"; //signal region preselection
  TString preselCutsZll = "(preselWord & (1<<2))"; //Z->ll control preselection  
  TString preselCutsWlv = "(preselWord & (1<<1))"; //W->lnu control preselection
  TString preselCutsPj = "(preselWord & (1<<5))"; //Photon+jets control preselection
  
  // monoJet selection
  TString metRawCut200 = "(metRaw > 200)";
  TString jetCuts = "((jet1.Pt() > 110 && abs(jet1.eta()) < 2.5))";
  TString dPhij1j2 = "TMath::ACos(TMath::Cos(jet2.phi() - jet1.phi()))";
  TString nJetCuts = TString("(njets == 1 || (njets == 2 && (")+dPhij1j2+TString(" < 2.0)))");
  TString tauVeto = "(ntaus == 0)"; //veto on taus
  TString photonVeto = "(nphotons == 0)"; //veto on photons 
  TString lepVeto = "(nlep == 0)"; //veto on leptons

  // Inclusive selection
  TString inclusiveCut = "(jet1.Pt() > 150 && abs(jet1.Eta()) < 2.0 && jet1CHF > 0.2 && jet1NHF < 0.7 && jet1NEMF < 0.7)";
  
  // BoostedV selection
  TString metRawCut250 = "(metRaw > 250)";
  TString fjetCut = "(fjet1.Pt() > 250 && abs(fjet1.Eta()) < 2.5)";
  TString bdt_loose = "(fjet1Tau2/fjet1Tau1 < 0.5 && 60 < fjet1MassPruned && fjet1MassPruned < 110)";
  //TString bdt_loose = "(bdt_noMass > 0.6)";
  //TString bdt_loose = "(bdt_pTOverM > 0.6)";
  TString bdt_tight = "(bdt_all > -0.1)";
  TString dPhifj1j2 = "TMath::ACos(TMath::Cos(jet2.phi() - fjet1.phi()))";
  TString fjetMassCut = "(fjet1.M() < 200.0)";
  TString sdCut = "(fjet1MassSDbm1 > 10)";
  TString prunedMassCut = ("(60 < fjet1MassPruned && fjet1MassPruned < 110)");
  TString qgCut = "(fjet1QGtag < 0.25)";
  TString nSubCut = "(fjet1Tau2/fjet1Tau1 < 0.5)";

  // Z->ll Control Region selection
  TString ZMuCut = "(genVdaughterId == 13)"; //require Z to decay to muons
  TString DiMuonCut     = "(nlep == 2 && (abs(lid1)>=13 && abs(lid2)>=13))"; //require two muons 
  TString InvMassCut = "(60 < mll && mll < 120)"; // require dilepton invariant mass to be within 30 GeV of Z mass
  TString metCorrectedZll = "TMath::Sqrt(TMath::Power(metRaw*TMath::Cos(metRawPhi) + lep1.Px() + lep2.Px(),2) + TMath::Power(metRaw*TMath::Sin(metRawPhi) + lep1.Py() + lep2.Py(),2))"; // MET calculated as if both muons were invisible
  TString metCorrectedZllCut200 = "(met > 200)"; // met corrected cut
  TString metCorrectedZllCut250 = "(met > 250)"; // met corrected cut
  
  // W->lv Control Region selection
  TString WMuCut = "(genVdaughterId == 13)"; // require W to decay to muon
  TString SingleMuonCut = "(nlep == 1 && abs(lid1) >= 13 && lep1.Pt() > 20 && abs(lep1.Eta()) < 2.4)"; // require one well defined muon
  TString JetLepdRWlv = "TMath::Sqrt(TMath::Power(fjet1.Eta() - lep1.Eta(),2)+TMath::Power(TMath::ACos(TMath::Cos(fjet1.phi() - lep1.Eta())),2))";
  TString JetLepdRWlvCut = TString("(")+JetLepdRWlv+TString(" > 0.3)");
  TString dPhiLepMet = "lep1.phi() - metRawPhi";
  TString TransMassCut ="(10 < mt && mt < 200)"; // require lepton tranverse mass to be W like
  TString metCorrectedWlvCut200 = "(met > 200)"; // met corrected cut
  TString metCorrectedWlvCut250 = "(met > 250)"; // met corrected cut

  // Photon+jets Control Region selection
  TString SinglePhotonCut = "(pho1.Pt() > 160 && abs(pho1.Eta()) < 2.5)"; // require one well defined photon
  TString metCorrectedPjCut200 = "(met > 200)"; // met corrected cut
  TString metCorrectedPjCut250 = "(met > 250)"; // met corrected cut
      
  // our plot task
  TString regions[4] = {"Met","Zll","Wlv","Pj"};
  TString cuts[14];
  
  TString MonoJetSelection           = eventWeight+monoJetTrigger+andC+jetCuts+andC+nJetCuts+andC+tauVeto+andC+photonVeto+andC+metFiltersCut;
  TString BoostedVVetoL              = "!("+fjetCut+andC+bdt_loose+")"; //BoostedV BDT based loose ANTI selection
  TString BoostedVSelectionL         = fjetCut+andC+bdt_loose; //BoostedV BDT based selection
  TString BoostedVSelectionT         = fjetCut+andC+bdt_tight; //BoostedV MVA tight selection
  TString SignalRegionSelection      = preselCutsZnunu+andC+lepVeto+andC+metRawCut200;
  TString ZllControlRegionSelection  = preselCutsZll+andC+DiMuonCut+andC+InvMassCut+andC+metCorrectedZllCut250;
  TString WlvControlRegionSelection = preselCutsWlv+andC+SingleMuonCut+andC+TransMassCut+andC+metCorrectedWlvCut250;
  TString PjControlRegionSelection = eventWeight+singlePhTrigger+andC+preselCutsPj+andC+lepVeto
                                    +andC+jetCuts+andC+nJetCuts+andC+tauVeto+andC+metFiltersCut
                                    +andC+SinglePhotonCut+andC+metCorrectedPjCut250;
  TString InclusiveZllControlRegionSelection  = preselCutsZll+andC+DiMuonCut+andC+InvMassCut+andC+metCorrectedZllCut200;
  TString InclusiveWlvControlRegionSelection = preselCutsWlv+andC+SingleMuonCut+andC+TransMassCut+andC+metCorrectedWlvCut200;
  TString InclusivePjControlRegionSelection = eventWeight+singlePhTrigger+andC+preselCutsPj+andC+lepVeto
                                    +andC+jetCuts+andC+nJetCuts+andC+tauVeto+andC+metFiltersCut
                                    +andC+SinglePhotonCut+andC+metCorrectedPjCut200;

  cuts[0]     = MonoJetSelection+andC+BoostedVSelectionL+andC+SignalRegionSelection+")";      // boostedV signal region selection, MET 200
  cuts[1]     = MonoJetSelection+andC+BoostedVSelectionT+andC+SignalRegionSelection+")";      // boostedV signal region selection, MET 200
  cuts[2]     = MonoJetSelection+andC+BoostedVSelectionL+andC+ZllControlRegionSelection+")";  // boostedV Z->ll control region selection, MET 250
  cuts[3]     = MonoJetSelection+andC+BoostedVSelectionL+andC+WlvControlRegionSelection+")";  // boostedV W->lnu control region selection, MET 250
  cuts[4]     = PjControlRegionSelection+andC+BoostedVSelectionL+")";                         // boostedV Photon+jets control region selection, MET 250
  cuts[5]     = MonoJetSelection+andC+BoostedVSelectionL+andC+metRawCut250+andC+SignalRegionSelection+")";        // boostedV signal region selection, MET 250
  cuts[6]     = MonoJetSelection+andC+inclusiveCut+andC+BoostedVVetoL+andC+SignalRegionSelection+")";             // inclusive signal region selection, MET 200
  cuts[7]     = MonoJetSelection+andC+inclusiveCut+andC+BoostedVVetoL+andC+InclusiveZllControlRegionSelection+")";// inclusive Z->ll region selection, MET 200
  cuts[8]     = MonoJetSelection+andC+inclusiveCut+andC+BoostedVVetoL+andC+InclusiveWlvControlRegionSelection+")";// inclusive W->lnu region selection, MET 200
  cuts[9]     = InclusivePjControlRegionSelection+andC+inclusiveCut+andC+BoostedVVetoL+")";                       // inclusive Photon+jets control region selection, MET 200
  cuts[10]    = MonoJetSelection+andC+inclusiveCut+andC+SignalRegionSelection+")";             // inclusive signal region selection, MET 200, no BDT cut
  cuts[11]    = MonoJetSelection+andC+inclusiveCut+andC+InclusiveZllControlRegionSelection+")";// inclusive Z->ll region selection, MET 200, no BDT cut
  cuts[12]    = MonoJetSelection+andC+inclusiveCut+andC+InclusiveWlvControlRegionSelection+")";// inclusive W->lnu region selection, MET 200, no BDT cut
  cuts[13]    = InclusivePjControlRegionSelection+andC+inclusiveCut+")";                       // inclusive Photon+jets control region selection, MET 200, no BDT cut 
  const int icut_no_BDT = 6;
  // Fix the region name if cut does not contain BDT
  if (icut >= icut_no_BDT)
    for (int i = 0; i < 4; i++)
      regions[i] = regions[i] + "_noVTagCut";
 
  PlotTask *plotTask = 0;
  

  // Prepare variable binning for met and jet pt
  double metbins150 [] = { 150.0,
                           160.0, 170.0, 180.0, 190.0, 200.0, 210.0, 220.0, 230.0, 240.0, 250.0,
                           260.0, 270.0, 280.0, 290.0, 300.0, 310.0, 320.0, 330.0, 350.0, 380.0,
                           430.0, 500.0, 1000.0 };
  double metbins200 [] = { 200.0, 210.0, 220.0, 230.0, 240.0, 250.0,
                           260.0, 270.0, 280.0, 290.0, 300.0, 310.0, 320.0, 330.0, 350.0, 380.0,
                           430.0, 500.0, 1000.0 };
  // Prepare variable binning for met and jet pt : boosted
  double metbins250[] = { 250.0, 300.0, 350.0, 400.0, 450.0, 500.0, 1000.};
  
  // Choose the correct binning 
  int nBinsVar = 0;
  if ( varbins && fabs(min-150) < 0.1 )
    nBinsVar =  sizeof(metbins150)/sizeof(double);
  if ( varbins && fabs(min-200) < 0.1 )
    nBinsVar =  sizeof(metbins200)/sizeof(double);
  if ( varbins && fabs(min-250) < 0.1 )
    nBinsVar =  sizeof(metbins250)/sizeof(double);

  double metbins[nBinsVar];
  if ( varbins && fabs(min-150) < 0.1 )
    for (int i=0; i< nBinsVar; i++) 
      metbins[i] =  metbins150[i];
  if ( varbins && fabs(min-200) < 0.1 )
    for (int i=0; i< nBinsVar; i++) 
      metbins[i] =  metbins200[i];
  if ( varbins && fabs(min-250) < 0.1 )
    for (int i=0; i< nBinsVar; i++) 
      metbins[i] =  metbins250[i];
    
  // For variable binning
  if (mode == -1) {
    gSystem->Setenv("MIT_ANA_CFG","boostedv-plots");
    // plot metRaw, BDT loose
    plotTask = new PlotTask(0,lumi);
    plotTask->SetVarBins(true);
    plotTask->SetHistRanges(metbins,min,max,0.,0.);
    int nBins = nBinsVar - 1;
    plotTask->SetNBins(nBins);
    plotTask->SetAxisTitles("E_{T}^{miss} [GeV]","Number of Events");
    plotTask->SetPngFileName("/tmp/dummy.png");
    plotTask->Plot(Stacked,nTuple,"metRaw",cuts[0],"BDT_loose_varbins_" + regions[0]);
    printf("Finished Plot for Region %s and Variable %s!\n",regions[0].Data(),"metRaw");
    delete plotTask;
  }
  // For limits
  if (mode == 0) {
    gSystem->Setenv("MIT_ANA_CFG","boostedv-plots");
    // plot metRaw, BDT loose
    plotTask = new PlotTask(0,lumi);
    plotTask->SetHistRanges(200.0,800.0,0.,0.);
    plotTask->SetNBins(60);
    plotTask->SetAxisTitles("E_{T}^{miss} [GeV]","Number of Events");
    plotTask->SetPngFileName("/tmp/dummy.png");
    plotTask->Plot(Stacked,nTuple,"metRaw",cuts[0],"BDT_loose_" + regions[0]);
    printf("Finished Plot for Region %s and Variable %s!\n",regions[0].Data(),"metRaw");
    delete plotTask;
  
    // plot metRaw, BDT tight
    plotTask = new PlotTask(0,lumi);
    plotTask->SetHistRanges(200.0,800.0,0.,0.);
    plotTask->SetNBins(60);
    plotTask->SetAxisTitles("E_{T}^{miss} [GeV]","Number of Events");
    plotTask->SetPngFileName("/tmp/dummy.png");
    plotTask->Plot(Stacked,nTuple,"metRaw",cuts[1],"BDT_tight_" + regions[0]);
    printf("Finished Plot for Region %s and Variable %s!\n",regions[0].Data(),"metRaw");
    delete plotTask;
  }
  // For signal region plots
  else if (mode == 1) {
    gSystem->Setenv("MIT_ANA_CFG","boostedv-plots");
    // plot metRaw
    plotTask = new PlotTask(0,lumi);
    plotTask->SetVarBins(varbins);
    plotTask->SetNormToWidth(varbins);
    if (!varbins)
      plotTask->SetHistRanges(min,max,0.,0.);
    else {
      plotTask->SetHistRanges(metbins,min,max,0.,0.);
      nBins = nBinsVar - 1;
    }
    plotTask->SetNBins(nBins);
    plotTask->SetAxisTitles(variable,"Number of Events");
    plotTask->SetPngFileName("/tmp/dummy.png");
    plotTask->Plot(Stacked,nTuple,variable,cuts[icut],"BDT_" + regions[0]);
    printf("Finished Plot for Region %s and Variable %s!\n",regions[0].Data(),variable.Data());
    delete plotTask;
  }
  // For Zll control region plots
  else if (mode == 2) {
    gSystem->Setenv("MIT_ANA_CFG","boostedv-plots-z");
    // plot met corrected
    plotTask = new PlotTask(0,lumi);
    plotTask->SetVarBins(varbins);
    plotTask->SetNormToWidth(varbins);
    if (!varbins)
      plotTask->SetHistRanges(min,max,0.,0.);
    else {
      plotTask->SetHistRanges(metbins,min,max,0.,0.);
      nBins = nBinsVar - 1;
    }
    plotTask->SetNBins(nBins);
    plotTask->SetAxisTitles(variable,"Number of Events");
    plotTask->SetPngFileName("/tmp/dummy.png");
    if (variable == "met")
      plotTask->Plot(Stacked,nTuple,metCorrectedZll,cuts[icut],"BDT_" + regions[1]);
    else 
      plotTask->Plot(Stacked,nTuple,variable,cuts[icut],"BDT_" + regions[1]);
    printf("Finished Plot for Region %s and Variable %s!\n",regions[1].Data(),variable.Data());
    delete plotTask;
  }
  // For Wlv control region plots
  else if (mode == 3) {
    gSystem->Setenv("MIT_ANA_CFG","boostedv-plots-w");
    // plot met corrected
    plotTask = new PlotTask(0,lumi);
    plotTask->SetVarBins(varbins);
    plotTask->SetNormToWidth(varbins);
    if (!varbins)
      plotTask->SetHistRanges(min,max,0.,0.);
    else {
      plotTask->SetHistRanges(metbins,min,max,0.,0.);
      nBins = nBinsVar - 1;
    }
    plotTask->SetNBins(nBins);
    plotTask->SetAxisTitles(variable,"Number of Events");
    plotTask->SetPngFileName("/tmp/dummy.png");
    if (variable == "met")
      plotTask->Plot(Stacked,nTuple,"met",cuts[icut],"BDT_" + regions[2]);
    else 
      plotTask->Plot(Stacked,nTuple,variable,cuts[icut],"BDT_" + regions[2]);
    printf("Finished Plot for Region %s and Variable %s!\n",regions[2].Data(),variable.Data());
    delete plotTask;
  }
  // For Photon+jets control region plots
  else if (mode == 4) {
    gSystem->Setenv("MIT_ANA_CFG","boostedv-plots-pj");
    // plot met corrected
    plotTask = new PlotTask(0,lumi);
    plotTask->SetVarBins(varbins);
    plotTask->SetNormToWidth(varbins);
    if (!varbins)
      plotTask->SetHistRanges(min,max,0.,0.);
    else {
      plotTask->SetHistRanges(metbins,min,max,0.,0.);
      nBins = nBinsVar - 1;
    }
    plotTask->SetNBins(nBins);
    plotTask->SetAxisTitles(variable,"Number of Events");
    plotTask->SetPngFileName("/tmp/dummy.png");
    if (variable == "met")
      plotTask->Plot(Stacked,nTuple,"met",cuts[icut],"BDT_" + regions[3]);
    else 
      plotTask->Plot(Stacked,nTuple,variable,cuts[icut],"BDT_" + regions[3]);
    printf("Finished Plot for Region %s and Variable %s!\n",regions[3].Data(),variable.Data());
    delete plotTask;
  }
  else 
    return;
  
  return;
}

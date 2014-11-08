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
void makeRawPlot(double lumi = 19700.0)
{
  // Create raw plots for further analysis (limits) or plotting

  // set the folder containing the input ntuples properly
  // here you can change the plot sources, i.e. the list of samples to be processed
  // and the location of the input flat trees
  gSystem->Setenv("MIT_ANA_CFG","boostedv-plots");
  gSystem->Setenv("MIT_ANA_HIST","/mnt/hscratch/dimatteo/boostedv-v8/merged/");

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
  TString preselCutsWlnu = "(preselWord & (1<<1))"; //W->lnu control preselection
  
  // monoJet selection
  TString metRawCut200 = "(metRaw > 200)";
  TString jetCuts = "((jet1.Pt() > 110 && abs(jet1.eta()) < 2.5))";
  TString dPhij1j2 = "TMath::ACos(TMath::Cos(jet2.phi() - jet1.phi()))";
  TString nJetCuts = TString("(njets == 1 || (njets == 2 && (")+dPhij1j2+TString(" < 2.5)))");
  TString tauVeto = "(ntaus == 0)"; //veto on taus
  TString photonVeto = "(nphotons == 0)"; //veto on photons 
  TString lepVeto = "(nlep == 0)"; //veto on leptons
  
  // BoostedV selection
  TString metRawCut250 = "(metRaw > 250)";
  TString fjetCut = "(fjet1.Pt() > 250 && abs(fjet1.Eta()) < 2.5)";
  TString bdt_all = "(bdt_all > -0.5)";
  TString dPhifj1j2 = "TMath::ACos(TMath::Cos(jet2.phi() - fjet1.phi()))";
  TString fjetMassCut = "(fjet1.M() < 200.0)";
  TString sdCut = "(fjet1MassSDbm1 > 10)";
  TString prunedMassCut = ("(65 < fjet1MassPruned && fjet1MassPruned < 105)");
  TString qgCut = "(fjet1QGtag < 0.25)";
  TString nSubCut = "(fjet1Tau2/fjet1Tau1 < 0.55)";

  // Z->ll Control Region selection
  TString ZMuCut = "(genVdaughterId == 13)"; //require Z to decay to muons
  TString DiMuonCut     = "(nlep == 2 && (lid1==13 && lid2==13))"; //require two muons 
  TString InvMass = "(TMath::Sqrt(TMath::Power(lep1.E() + lep2.E(), 2) - TMath::Power(lep1.Px() + lep2.Px(), 2) - TMath::Power(lep1.Py() + lep2.Py(), 2) - TMath::Power(lep1.Pz() + lep2.Pz(), 2)))"; // dilepton invariant mass
  TString InvMassCut = TString("(60 < ")+InvMass+TString(") && (")+InvMass+TString(" < 120)"); // require dilepton invariant mass to be within 30 GeV of Z mass
  TString metCorrectedZll = "TMath::Sqrt(TMath::Power(metRaw*TMath::Cos(metRawPhi) + lep1.Px() + lep2.Px(),2) + TMath::Power(metRaw*TMath::Sin(metRawPhi) + lep1.Py() + lep2.Py(),2))"; // MET calculated as if both muons were invisible
  TString metCorrectedZllCut200 = TString("(")+metCorrectedZll+TString(" > 200)"); // met corrected cut
  TString metCorrectedZllCut250 = TString("(")+metCorrectedZll+TString(" > 250)"); // met corrected cut

  
  // W->lv Control Region selection
  TString WMuCut = "(genVdaughterId == 13)"; // require W to decay to muon
  TString SingleMuonCut = "(nlep == 1 && abs(lid1) == 13 && lep1.Pt() > 20 && abs(lep1.Eta()) < 2.4)"; // require one well defined muon
  TString JetLepdRWlnu = "TMath::Sqrt(TMath::Power(fjet1.Eta() - lep1.Eta(),2)+TMath::Power(TMath::ACos(TMath::Cos(fjet1.phi() - lep1.Eta())),2))";
  TString JetLepdRWlnuCut = TString("(")+JetLepdRWlnu+TString(" > 0.3)");
  TString dPhiLepMet = "lep1.phi() - metRawPhi";
  TString TransMass = TString("(TMath::Sqrt(2*metRaw*lep1.Pt()*(1-TMath::Cos(")+dPhiLepMet+TString("))))"); // lepton tranverse mass
  TString TransMassCut = TString("(10 < ")+TransMass+TString(") && (")+TransMass+TString(" < 200)"); // require lepton tranverse mass to be W like
  TString metCorrectedWlnu = "TMath::Sqrt(TMath::Power(metRaw*TMath::Cos(metRawPhi) + lep1.Px(),2) + TMath::Power(metRaw*TMath::Sin(metRawPhi) + lep1.Py(),2))"; // MET calculated as if both muons were invisible
  TString metCorrectedWlnuCut200 = TString("(")+metCorrectedWlnu+TString(" > 200)"); // met corrected cut
  TString metCorrectedWlnuCut250 = TString("(")+metCorrectedWlnu+TString(" > 250)"); // met corrected cut
      
  // our plot task
  TString regions[6] = {"Monojet_Signal","Monojet_Zll","Monojet_Wlnu","BoostedV_Signal","BoostedV_Zll","BoostedV_Wlnu"};
  TString cuts[6];
  
  TString MonoJetSelection           = eventWeight+monoJetTrigger+andC+nJetCuts+andC+tauVeto+andC+photonVeto+andC+metFiltersCut;
  TString BoostedVSelection          = fjetCut+andC+bdt_all; //BoostedV MVA based selection
  TString SignalRegionSelection      = preselCutsZnunu+andC+lepVeto+andC+metRawCut200;
  TString ZllControlRegionSelection  = preselCutsZll+andC+DiMuonCut+andC+InvMassCut+andC+metCorrectedZllCut200;
  TString WlnuControlRegionSelection = preselCutsWlnu+andC+SingleMuonCut+andC+TransMassCut+andC+metCorrectedWlnuCut200;

  cuts[0]     = MonoJetSelection+andC+SignalRegionSelection+")";       // monojet signal region selection
  cuts[1]     = MonoJetSelection+andC+ZllControlRegionSelection+")";   // monojet Z->ll control region selection
  cuts[2]     = MonoJetSelection+andC+WlnuControlRegionSelection+")";  // monojet W->lnu control region selection
  cuts[3]     = MonoJetSelection+andC+BoostedVSelection+andC+SignalRegionSelection+")";      // boostedV signal region selection
  cuts[4]     = MonoJetSelection+andC+BoostedVSelection+andC+ZllControlRegionSelection+")";  // boostedV Z->ll control region selection
  cuts[5]     = MonoJetSelection+andC+BoostedVSelection+andC+WlnuControlRegionSelection+")"; // boostedV W->lnu control region selection

 
  // plot metRaw
  PlotTask *plotTask = new PlotTask(0,lumi);
  plotTask->SetHistRanges(200.0,800.0,0.,0.);
  plotTask->SetNBins(60);
  plotTask->SetAxisTitles("E_{T}^{miss} [GeV]","Number of Events");
  plotTask->SetPngFileName("/tmp/dummy.png");
  plotTask->Plot(Stacked,nTuple,"metRaw",cuts[3],"BDT");
  printf("Finished Plot for Region %s and Variable %s!\n",regions[3].Data(),"metRaw");
  
  return;
}

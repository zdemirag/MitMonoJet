//--------------------------------------------------------------------------------------------------
// Reduce the ntuples and compute all the necessary weights for further studies and plots
//
// Authors: L. Di Matteo                                                                  (Dec 2014)
//--------------------------------------------------------------------------------------------------
#include <iostream>
#include <TSystem.h>
#include <TROOT.h>
#include <TFile.h>
#include <TNtuple.h>
#include <TTree.h>
#include <TH1D.h>
#include <TMath.h>
#include "MitCommon/MathTools/interface/MathUtils.h"
#include "MitPlots/Input/interface/TaskSamples.h"
#include "MitMonoJet/Core/MitDMSTree.h"
#include "MitMonoJet/Core/MitLimitTree.h"

#include "TLorentzVector.h"

using namespace std;
using namespace mithep;

//---
TString getEnv(const char* name);
//---
void fillOutNtuples(MitLimitTree &outtree, MitDMSTree &intree, double baseWeight, int selMode = 0, bool isData = false, bool exclusive = true, bool testing = false);
//---
bool eventPassSelection(MitDMSTree &intree, int selMode = 0, bool exclusive = true);
//==================================================================================================
void makeReducedTree(int selMode = 0, double lumi = 19700.0, bool updateFile = false, bool exclusive = true, bool testing = false) 
{
  // Define tree name (depends on selection)
  TString outTreeNameExt;
  if (selMode == 0 || selMode == 4 || selMode == 8) 
    outTreeNameExt = "_signal";
  else if (selMode == 1 || selMode == 5 || selMode == 9) 
    outTreeNameExt = "_di_muon_control";
  else if (selMode == 2 || selMode == 6 || selMode == 10) 
    outTreeNameExt = "_single_muon_control";
  else if (selMode == 3 || selMode == 7 || selMode == 11) 
    outTreeNameExt = "_photon_control";
  else 
    cout << "ERROR -- Incorrect selMode parameter, please review!" << endl;

  // Define output files mode
  TString outFileMode = "RECREATE";
  if (updateFile)
    outFileMode = "UPDATE";

  // Read all environment variables
  TString anaDir = getEnv("MIT_MONOJET_DIR");
  //TString hstDir = getEnv("MIT_ANA_HIST");
  TString hstDir = "/scratch4/dimatteo/cms/hist/boostedv-v10/merged-test/";  
  TString anaCfg = "boostedv-ana";
  //TString prdCfg = getEnv("MIT_PROD_CFG");
  TString prdCfg = "boostedv-v10";
  
  // Fix data list for photons
  if (selMode == 3 || selMode == 7 || selMode == 11)
    anaCfg = "boostedv-ana-pj";

  if (testing)
    anaCfg = "boostedv-ana-test";

  // Define samples
  TaskSamples* samples = new TaskSamples(prdCfg.Data(),hstDir.Data());
  samples->SetNameTxt(anaCfg.Data());
  samples->ReadFile((anaDir + TString("/config")).Data());
  vector<const Sample*> listOfSamples;
  for (UInt_t iSample=0; iSample < samples->NDataSamples(); iSample++) listOfSamples.push_back(samples->GetDataSample(iSample));
  for (UInt_t iSample=0; iSample < samples->NSamples(); iSample++) listOfSamples.push_back(samples->GetSample(iSample));  

  // Prepare pointer to outfile
  TFile *fin;
  TFile *fout;
  TString outFileName = "boosted.root";
  if (selMode >= 4 && selMode < 8)
    outFileName = "resolved.root";
  if (selMode >= 8)
    outFileName = "inclusive.root";
  if (!exclusive)
    outFileName = "baseline.root";
  if (testing)
    outFileName = "testing.root";
  fout = new TFile(outFileName,outFileMode);
  
  // Prepare object to store outtree
  MitLimitTree outTree;
  
  // Generate reduced trees
  // loop through the samples and produce the reduced trees
  for (UInt_t iSample=0; iSample < listOfSamples.size(); iSample++) {
    TString inFilePath = hstDir+"/"+*listOfSamples.at(iSample)->File();
    fin = new TFile(inFilePath,"READ");
    cout << "INFO -- Reading sample file: " << *listOfSamples.at(iSample)->File() << endl;
    // Prepare event weight
    double thisXsec = *listOfSamples.at(iSample)->Xsec();
    double nGenEvts = ((TH1D*)fin->FindObjectAny("hDAllEvents"))->GetEntries();
    double baseWeight = lumi*thisXsec/nGenEvts;
    bool isData = false;
    if (baseWeight < 0)
      isData = true;
    // Read input tree
    MitDMSTree inTree;
    inTree.LoadTree(inFilePath,1);
    inTree.InitTree(1);
    // Start a new sample group according to cfg file legend
    if (*listOfSamples.at(iSample)->Legend() != " ") {
      // Close previous group
      if (iSample > 0) {
        fout->cd();
        outTree.tree_->Write();
        fout->Close();
        fout = TFile::Open(outFileName,"UPDATE");     
      } 
      fout->cd();
      TString outTreeName = *listOfSamples.at(iSample)->Legend();
      outTreeName += outTreeNameExt;
      outTree.CreateTree(outTreeName.Data());
      outTree.InitVariables();
    }
    // Scan on input and fill output ntuple
    cout << "INFO ---> Number of events passing the selection is: ";
    fillOutNtuples(outTree,inTree,baseWeight,selMode,isData,exclusive,testing);
    
    // Close last group
    if (iSample == (listOfSamples.size()-1)) {
      fout->cd();
      outTree.tree_->Write(); 
      fout->Close();
      break;
    }      

    fin->Close();
  }  
  
  return;
}


//==================================================================================================
TString getEnv(const char* name)
{
  if (! gSystem->Getenv(name)) {
    printf(" Environment variable: %s  not defined. EXIT!\n",name);
    return TString("");
  } 
  return TString(gSystem->Getenv(name));  
}

//==================================================================================================
void fillOutNtuples(MitLimitTree &outtree, MitDMSTree &intree, double baseWeight, int selMode, bool isData, bool exclusive, bool testing)
{
  double weight = -1;
  double sumweight = 0;
  // Loop over tree entries
  Int_t nEntries = intree.tree_->GetEntries();
  for ( Int_t iEntry = 0; iEntry < nEntries; iEntry++ ) {
    
    // Get this tree entry
    intree.tree_-> GetEntry(iEntry);

    // Determine if event passes selection
    if (!eventPassSelection(intree,selMode,exclusive))
      continue;
      
    // Determine correctly the event weights
    if (!isData)
      weight = baseWeight*intree.puweight_;    
    sumweight += weight;
    
    // Determine the outtree variables for this event
    outtree.mvamet_ = intree.met_;
    outtree.mvametphi_ = intree.metPhi_;
    outtree.jet1pt_ = intree.fjet1_.Pt();
    outtree.genjetpt_ = intree.fjet1_.Pt();
    if (selMode >= 4 && selMode < 8)
      outtree.jet1pt_ = (intree.rjet1_ + intree.rjet2_).Pt();
      outtree.genjetpt_ = intree.rjet1_.Pt();
    if (selMode >= 8)
      outtree.jet1pt_ = intree.jet1_.Pt();
      outtree.genjetpt_ = intree.jet1_.Pt();
    outtree.genVpt_ = intree.genV_.Pt();
    outtree.weight_ = weight;
    if (isData) {
      outtree.genjetpt_ = -1.;
      outtree.genVpt_ = -1.;
      outtree.weight_ = -1.;
    }

    // Fill output tree
    outtree.tree_->Fill();
  }
  
  if (isData)
    sumweight = -sumweight;
  if (!isData || testing)
    cout << sumweight << endl;
  else 
    cout << "BLINDED!" << endl;
  return;
}

//==================================================================================================
bool eventPassSelection(MitDMSTree &intree, int selMode, bool exclusive)
{
  // Trigger
  bool triggerBit = ((intree.trigger_ & (1<<0)) || (intree.trigger_ & (1<<1)));
  if (selMode == 3 || selMode == 7 || selMode == 11)
    triggerBit = (intree.trigger_ & (1<<3));

  // Met filters
  bool metFiltersBit = (intree.metFiltersWord_ == 511 || intree.metFiltersWord_ == 1023);
  
  // Preselection
  bool preselBit = (intree.preselWord_ & (1<<3)); //signal region preselection
  if (selMode == 1 || selMode == 5 || selMode == 9)
    preselBit = (intree.preselWord_ & (1<<2)); //Z->ll control preselection
  if (selMode == 2 || selMode == 6 || selMode == 10)
    preselBit = (intree.preselWord_ & (1<<1)); //W->lnu control preselection
  if (selMode == 3 || selMode == 7 || selMode == 11)
    preselBit = (intree.preselWord_ & (1<<5)); //Photon+jets control preselection

  // Met
  bool metBit = (intree.met_ > 180. && intree.met_ < 1000.);
  if (selMode < 4)
    metBit = (intree.met_ > 200.);
                      
  // Narrow jets
  bool jetBit = ((intree.jet1_.Pt() > 110 && abs(intree.jet1_.eta()) < 2.5)
               && intree.jet1CHF_ > 0.2 && intree.jet1NHF_ < 0.7 && intree.jet1NEMF_ < 0.7);
  bool nJetBit = (intree.njets_ == 1 || (intree.njets_ == 2 && 
                  abs(MathUtils::DeltaPhi(intree.jet1_.phi(),intree.jet2_.phi())) < 2.0));

  // Resolved category
  bool resolvedBit = (intree.rmvaval_ > 0.6);

  // Inclusive category
  bool inclusiveBit = (intree.jet1_.Pt() > 150 && abs(intree.jet1_.Eta()) < 2.0);

  // Fat jet category::cut-based
  bool fatJetBit = (intree.fjet1_.Pt() > 200 && abs(intree.fjet1_.Eta()) < 2.5);
  fatJetBit = fatJetBit && (intree.fjet1Tau2_/intree.fjet1Tau1_  < 0.5 
                            && 60 < intree.fjet1MassPruned_ && intree.fjet1MassPruned_ < 110
                            && intree.met_ > 200.);

  // Vetoes
  bool vetoBit = (intree.ntaus_ == 0);
  if (selMode == 0 || selMode == 4 || selMode == 8)
    vetoBit = vetoBit && (intree.nphotons_ == 0 && intree.nlep_ == 0);
  if (selMode == 1 || selMode == 5 || selMode == 9)
    vetoBit = vetoBit && (intree.nphotons_ == 0);
  if (selMode == 2 || selMode == 6 || selMode == 10)
    vetoBit = vetoBit && (intree.nphotons_ == 0);
  if (selMode == 3 || selMode == 7 || selMode == 11)
    vetoBit = vetoBit && (intree.nlep_ == 0);
    
  // Extra
  bool extraBit = true;
  if (selMode == 1 || selMode == 5 || selMode == 9) {
    extraBit = extraBit && (intree.nlep_ == 2 && (abs(intree.lid1_)>=130 && abs(intree.lid2_)>=13)
                                              && (abs(intree.lep1_.Eta()) < 2.1 && abs(intree.lep2_.Eta()) < 2.1));
    extraBit = extraBit && (60 < intree.mll_ && intree.mll_ < 120);    
  } //Zll
  if (selMode == 2 || selMode == 6 || selMode == 10) {
    extraBit = extraBit && 
              (intree.nlep_ == 1 && abs(intree.lid1_) >= 1300 && intree.lep1_.Pt() > 20 && abs(intree.lep1_.Eta()) < 2.1);
    extraBit = extraBit && (intree.mt_ > 40. && intree.mt_ < 200.);
  } //Wlv
  if (selMode == 3 || selMode == 7 || selMode == 11) {
    extraBit = extraBit && 
              (intree.pho1_.Pt() > 160 && abs(intree.pho1_.Eta()) < 2.5);
  } //Pj
        
  //cout 
  //<< triggerBit << " " << metFiltersBit << " " << preselBit << " " << metBit << " " 
  //<< jetBit << " " << resolvedBit << " " << inclusiveBit << " " <<  fatJetBit << " " << vetoBit << " " << extraBit 
  //<< endl;

  bool theDecision;
  // Boosted selection
  if (selMode < 4) 
    theDecision = triggerBit && metFiltersBit && preselBit && metBit && jetBit && nJetBit && fatJetBit && vetoBit && extraBit;
  // Resolved selection and discard boosted events
  else if (selMode >= 4 && selMode < 8)
    theDecision = triggerBit && metFiltersBit && preselBit && metBit && resolvedBit && vetoBit && extraBit && !fatJetBit;
  // Inclusive selection and discard boosted/resolved events
  else if (selMode >= 8 && exclusive)
    theDecision = triggerBit && metFiltersBit && preselBit && metBit && jetBit && nJetBit && inclusiveBit && vetoBit && extraBit && !fatJetBit && !resolvedBit;
  // Inclusive selection keep boosted/resolved events
  else if (selMode >= 4 && !exclusive)
    theDecision = triggerBit && metFiltersBit && preselBit && metBit && jetBit && nJetBit && inclusiveBit && vetoBit && extraBit;
  else 
    cout << "ERROR - Selection logic is wrong! Please fix" << endl;
  
  return theDecision;
}

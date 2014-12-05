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

#include "TLorentzVector.h"

using namespace std;
using namespace mithep;

//---
TString getEnv(const char* name);
//---
void fillOutNtuples(TNtuple* ntuple, MitDMSTree &intree, double baseWeight, int selMode = 0, bool isData = false);
//---
bool eventPassSelection(MitDMSTree &intree, float &met, int selMode = 0);
//==================================================================================================
void makeReducedTree(int selMode = 0, double lumi = 19700.0, bool updateFile = false) 
{
  // Define tree name (depends on selection)
  TString ntupleName = "Process_signal";
  if (selMode == 1) 
    ntupleName = "Process_di_muon_control";
  else if (selMode == 2) 
    ntupleName = "Process_single_muon_control";
  else if (selMode == 3) 
    ntupleName = "Process_photon_control";
  else 
    cout << "ERROR -- Incorrect selMode parameter, please review!" << endl;

  // Define output files mode
  TString outFileMode = "RECREATE";
  if (updateFile)
    outFileMode = "UPDATE";

  // Read all environment variables
  TString anaDir = getEnv("MIT_MONOJET_DIR");
  TString hstDir = getEnv("MIT_ANA_HIST");
  TString anaCfg = "boostedv-ana";
  TString prdCfg = getEnv("MIT_PROD_CFG");
  
  // Fix data list for photons
  if (selMode == 3)
    anaCfg = "boostedv-ana-pj";

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
  
  // Prepare pointer to outntuple
  TNtuple *ntuple;
  
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
        ntuple->Write();
        fout->Close();
      }      
      fout = new TFile(*listOfSamples.at(iSample)->Legend()+".root",outFileMode);
      fout->cd();
      ntuple = new TNtuple(ntupleName,"Limit ntuple","mvamet:jet1pt:genjetpt:genVpt:weight");      
    }
    // Scan on input and fill output ntuple
    cout << "INFO ---> Number of events passing the selection is: ";
    fillOutNtuples(ntuple,inTree,baseWeight,selMode,isData);
    
    // Close last group
    if (iSample == (listOfSamples.size()-1)) {
    //if (iSample == 6) {
      fout->cd();
      ntuple->Write();
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
void fillOutNtuples(TNtuple* ntuple, MitDMSTree &intree, double baseWeight, int selMode, bool isData)
{
  double weight = -1;
  float met = -1;
  double sumweight = 0;
  // Loop over tree entries
  Int_t nEntries = intree.tree_->GetEntries();
  for ( Int_t iEntry = 0; iEntry < nEntries; iEntry++ ) {
    
    // Get this tree entry
    intree.tree_-> GetEntry(iEntry);

    // Determine if event passes selection
    if (!eventPassSelection(intree,met,selMode))
      continue;
      
    // Determine correctly the event weights
    if (!isData)
      weight = baseWeight*intree.puweight_;    
    sumweight += weight;

    // Fill output tree
    if (!isData)
      ntuple->Fill(met,intree.fjet1_.Pt(),intree.fjet1_.Pt(),intree.genV_.Pt(),weight);
    else 
      ntuple->Fill(met,intree.fjet1_.Pt(),-1.,-1.,-1.);
  }
  
  if (!isData)
    cout << sumweight << endl;
  else 
    cout << "BLINDED!" << endl;
  return;
}

//==================================================================================================
bool eventPassSelection(MitDMSTree &intree, float &met, int selMode)
{
  // Trigger
  bool triggerBit = ((intree.trigger_ & (1<<0)) || (intree.trigger_ & (1<<1)));
  if (selMode == 3)
    triggerBit = (intree.trigger_ & (1<<3));

  // Met filters
  bool metFiltersBit = (intree.metFiltersWord_ == 511 || intree.metFiltersWord_ == 1023);
  
  // Preselection
  bool preselBit = (intree.preselWord_ & (1<<3)); //signal region preselection
  if (selMode == 1)
    preselBit = (intree.preselWord_ & (1<<2)); //Z->ll control preselection
  if (selMode == 2)
    preselBit = (intree.preselWord_ & (1<<1)); //W->lnu control preselection
  if (selMode == 3)
    preselBit = (intree.preselWord_ & (1<<5)); //Photon+jets control preselection

  // Met
  met = intree.metRaw_;
  if (selMode == 1)
    met = TMath::Sqrt(TMath::Power(intree.metRaw_*TMath::Cos(intree.metRawPhi_) + 
                      intree.lep1_.Px() + intree.lep2_.Px(),2) + 
                      TMath::Power(intree.metRaw_*TMath::Sin(intree.metRawPhi_) + 
                      intree.lep1_.Py() + intree.lep2_.Py(),2)); //Z->ll control preselection
  if (selMode == 2)
    met = TMath::Sqrt(TMath::Power(intree.metRaw_*TMath::Cos(intree.metRawPhi_) + 
                      intree.lep1_.Px(),2) + 
                      TMath::Power(intree.metRaw_*TMath::Sin(intree.metRawPhi_) + 
                      intree.lep1_.Py(),2)); //W->lnu control preselection
  if (selMode == 3)
    met = TMath::Sqrt(TMath::Power(intree.metRaw_*TMath::Cos(intree.metRawPhi_) + 
                      intree.pho1_.Px(),2) + 
                      TMath::Power(intree.metRaw_*TMath::Sin(intree.metRawPhi_) + 
                      intree.pho1_.Py(),2)); //Photon+jets control preselection
  bool metBit = (met > 250.);
                      
  // Narrow jets
  bool jetBit = ((intree.jet1_.Pt() > 110 && abs(intree.jet1_.eta()) < 2.5));
  jetBit = jetBit && (intree.njets_ == 1 || (intree.njets_ == 2 && 
                      MathUtils::DeltaPhi(intree.jet1_.phi(),intree.jet2_.phi()) < 2.5));

  // Fat jet :: FIXME with bdt
  bool fatJetBit = (intree.fjet1_.Pt() > 250 && abs(intree.fjet1_.Eta()) < 2.5);
  fatJetBit = fatJetBit && (intree.bdt_all_ > -0.5);

  // Vetoes
  bool vetoBit = (intree.ntaus_ == 0);
  if (selMode == 0)
    vetoBit = vetoBit && (intree.nphotons_ == 0 && intree.nlep_ == 0);
  if (selMode == 1)
    vetoBit = vetoBit && (intree.nphotons_ == 0);
  if (selMode == 2)
    vetoBit = vetoBit && (intree.nphotons_ == 0);
  if (selMode == 3)
    vetoBit = vetoBit && (intree.nlep_ == 0);
    
  // Extra
  bool extraBit = true;
  if (selMode == 1) {
    extraBit = extraBit && (intree.nlep_ == 2 && (intree.lid1_==13 && intree.lid2_==13));
    TLorentzVector tempBos;
    TLorentzVector tempLep1;
    tempLep1.SetPtEtaPhiE(intree.lep1_.Pt(),intree.lep1_.Eta(),intree.lep1_.Phi(),intree.lep1_.E());
    TLorentzVector tempLep2;
    tempLep2.SetPtEtaPhiE(intree.lep2_.Pt(),intree.lep2_.Eta(),intree.lep2_.Phi(),intree.lep2_.E());
    tempBos = tempLep1 + tempLep2;
    extraBit = extraBit && (tempBos.M() > 60. && tempBos.M() < 120.);    
  } //Zll
  if (selMode == 2) {
    extraBit = extraBit && 
              (intree.nlep_ == 1 && abs(intree.lid1_) == 13 && intree.lep1_.Pt() > 20 && abs(intree.lep1_.Eta()) < 2.4);
    float mt = sqrt(2*intree.metRaw_*intree.lep1_.Pt()*(1-TMath::Cos(intree.metRawPhi_-intree.lep1_.Phi())));
    extraBit = extraBit && (mt > 10. && mt < 200.);
  } //Wlv
  if (selMode == 3) {
    extraBit = extraBit && 
              (intree.nphotons_ == 1 && intree.pho1_.Pt() > 160 && abs(intree.pho1_.Eta()) < 2.5);
  } //Pj
        
  //cout << triggerBit << " " << metFiltersBit << " " << preselBit << " " << metBit << " " << jetBit << " " << fatJetBit << " " << vetoBit << " " << extraBit << endl;
  bool theDecision = triggerBit && metFiltersBit && preselBit && metBit && jetBit && fatJetBit && vetoBit && extraBit;
  
  return theDecision;
}

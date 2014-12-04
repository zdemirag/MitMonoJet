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
#include "MitPlots/Input/interface/TaskSamples.h"
#include "MitMonoJet/Core/MitDMSTree.h"

using namespace std;
using namespace mithep;

//---
TString getEnv(const char* name);
//---
void fillOutNtuples(TNtuple* ntuple, MitDMSTree &intree, double baseWeight);
//==================================================================================================
void makeReducedTree(double lumi = 19700.0) 
{
  // Define tree name (depends on selection)
  TString ntupleName = "Process_signal";

  // Read all environment variables
  TString anaDir = getEnv("MIT_MONOJET_DIR");
  TString hstDir = getEnv("MIT_ANA_HIST");
  TString anaCfg = "boostedv-ana";
  TString prdCfg = getEnv("MIT_PROD_CFG");

  // Define samples
  TaskSamples* samples = new TaskSamples(prdCfg.Data(),hstDir.Data());
  samples->SetNameTxt(anaCfg.Data());
  samples->ReadFile((anaDir + TString("/config")).Data());
  vector<const Sample*> listOfSamples;
  //for (UInt_t iSample=0; iSample < samples->NDataSamples(); iSample++) listOfSamples.push_back(samples->GetDataSample(iSample));
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
    cout << "read in file" << endl;
    // Prepare event weight
    double thisXsec = *listOfSamples.at(iSample)->Xsec();
    double nGenEvts = ((TH1D*)fin->FindObjectAny("hDAllEvents"))->GetEntries();
    double baseWeight = lumi*thisXsec/nGenEvts;
    cout << "got weight" << endl;
    // Read input tree
    MitDMSTree inTree;
    inTree.LoadTree(inFilePath,0);
    inTree.InitTree(0);
    cout << "setup in tree" << endl;
    // Start a new sample group according to cfg file legend
    if (*listOfSamples.at(iSample)->Legend() != "~") {
      // Close previous group
      if (iSample > 0) {
        fout->cd();
        ntuple->Write();
        fout->Close();
      }      
      fout = new TFile(*listOfSamples.at(iSample)->Legend()+".root","RECREATE");
      fout->cd();
      ntuple = new TNtuple(ntupleName,"Limit ntuple","mvamet:jet1pt:genjetpt:genVpt:weight");      
      cout << "opened new ntuple" << endl;
    }
    // Scan on input and fill output ntuple
    fillOutNtuples(ntuple,inTree,baseWeight);
    cout << "called filler function" << endl;
    
    // Close last group
    //if (iSample == (listOfSamples.size()-1)) {
    if (iSample == 0) {
      fout->cd();
      ntuple->Write();
      fout->Close();
    }      

    fin->Close();
    break;
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
void fillOutNtuples(TNtuple* ntuple, MitDMSTree &intree, double baseWeight)
{
  double weight = -1;
  // Loop over tree entries
  Int_t nEntries = intree.tree_->GetEntries();
  for ( Int_t iEntry = 0; iEntry < nEntries; iEntry++ ) {
    
    // Get this tree entry
    intree.tree_-> GetEntry(iEntry);

    weight = baseWeight*intree.puweight_;

    // Fill output tree
    ntuple->Fill(intree.mvamet_,intree.fjet1_.Pt(),intree.fjet1_.Pt(),intree.genV_.Pt(),weight);
  }
  
  return;
}

#include <iostream>
#include <string>
#include "TFile.h"
#include "TTree.h"

void patchTree(std::string iFile0="../samples/boostedv-v8_s12-zll-ptz100-v7c_noskim_flatntuple.root",
               std::string iFile1="Output.root",std::string iName0="all") { 
  TFile *lFile0 = new TFile(iFile0.c_str(),"UPDATE");
  TTree *lTree0 = (TTree*) lFile0->FindObjectAny("DMSTree");

  TFile *lFile1 = new TFile(iFile1.c_str());
  TTree *lTree1 = (TTree*) lFile1->FindObjectAny("DMSTree");
  
  //if(iName0 == iName1 && iName1 != "2.*jet1qg2+jet1qg1") iName1 = "";
  if(iName0 == "jet1tau2/jet1tau1")  iName0 = "t2t1";
  if(iName0 == "2.*jet1qg2+jet1qg1") iName0 = "qg12";
  std::string iName = "bdt_"+iName0;
  // cout << " Test : " << iName0 << endl;

  float lBDT0 = 0; lTree1->SetBranchAddress(iName.c_str(),&lBDT0); 
  TBranch* lBranch = lTree1->GetBranch(iName.c_str());
  
  iName = "bdt_"+iName0;
  float lBDT  = 0; 
  TBranch* lBDTBranch = lTree0->Branch(iName.c_str(),&lBDT,(iName+"/F").c_str());
  //lTree0->Branch(iName.c_str(),&lBDT);
  //TBranch* lBDTBranch = lTree0->GetBranch(iName.c_str());
  int lNEntries = lTree0->GetEntriesFast();

  for(int i0 = 0; i0 < lNEntries; i0++) { 
    lBranch->GetEntry(i0); 
    lBDT = lBDT0;
    //lTree0->GetEntry(i0); 
    lBDTBranch->Fill();
  }
  lFile0->cd();
  lTree0->Write();
}

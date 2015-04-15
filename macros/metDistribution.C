#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TSystem.h>
#include "MitAna/TreeMod/interface/Analysis.h"
#include "MitMonoJet/Mods/interface/MetAnalysisMod.h"
#include "MitAna/DataTree/interface/PFJetCol.h"
#include "MitAna/PhysicsMod/interface/PublisherMod.h"
#include "MitPhysics/Mods/interface/GoodPVFilterMod.h"
#include "MitPhysics/Mods/interface/JetCorrectionMod.h"
#include "MitPhysics/Mods/interface/MetCorrectionMod.h"
#endif

void metDistribution(const char *fileset    = "0000",
             const char *skim       = "noskim",
             const char *dataset    = "s12-zjets-ptz100-v7a",
             const char *book       = "",
             const char * = "",
             const char *outputLabel = "monojet",
             int         nEvents    = 1000)
{
  TString jsonEnv(gSystem->Getenv("MIT_PROD_JSON"));
  bool isData = jsonEnv.Length() != 0 && jsonEnv != "~";

  TString MitData(gSystem->Getenv("MIT_DATA"));

  mithep::Analysis* ana = new mithep::Analysis;
  ana->SetUseHLT(kFALSE);
  if(nEvents >= 0)
    ana->SetProcessNEvents(nEvents);

  // File names for input and output are identical
  // so that default condor submission script works.
  // Should change the output file name to be transferred
  // in MitAna/bin/submit.sh if you use a different name.
  TString fileName(outputLabel);
  fileName += "_";
  fileName += dataset;
  fileName += "_";
  fileName += skim;
  fileName += "_";
  fileName += fileset;
  fileName += ".root";

  TString outputDirName;
  TString inputDirName;
  if(strlen(book) != 0){
    // Assume skim made by batch submission as input
    outputDirName = ".";

    UserGroup_t* userInfo = gSystem->GetUserInfo();
    inputDirName = "/mnt/hscratch/" + userInfo->fUser;
    inputDirName += "/skim/";
    inputDirName += outputLabel;
    inputDirName += "/";
    inputDirName += book;
    inputDirName += "/";
    inputDirName += dataset;
    if(TString(skim) != "noskim"){
      inputDirName += "/";
      inputDirName += skim;
    }

    delete userInfo;
  }
  else{
    outputDirName = gSystem->Getenv("MIT_PROD_HIST");    

    inputDirName += ".";
  }

  ana->SetOutputName(outputDirName + "/" + fileName);

  ana->AddFile(inputDirName + "/" + fileName);

  mithep::PublisherMod<mithep::PFJet, mithep::Jet> *pubJet =
    new mithep::PublisherMod<mithep::PFJet, mithep::Jet>("JetPub");
  pubJet->SetInputName("AKt5PFJets");
  pubJet->SetOutputName("PubAKt5PFJets");

  mithep::JetCorrectionMod *jetCorr = new mithep::JetCorrectionMod;
  if(isData){ 
    jetCorr->AddCorrectionFromFile((MitData+TString("/Summer13_V1_DATA_L1FastJet_AK5PF.txt")).Data()); 
    jetCorr->AddCorrectionFromFile((MitData+TString("/Summer13_V1_DATA_L2Relative_AK5PF.txt")).Data()); 
    jetCorr->AddCorrectionFromFile((MitData+TString("/Summer13_V1_DATA_L3Absolute_AK5PF.txt")).Data()); 
    jetCorr->AddCorrectionFromFile((MitData+TString("/Summer13_V1_DATA_L2L3Residual_AK5PF.txt")).Data());
  }                                                                                      
  else {                                                                                 
    jetCorr->AddCorrectionFromFile((MitData+TString("/Summer13_V1_MC_L1FastJet_AK5PF.txt")).Data()); 
    jetCorr->AddCorrectionFromFile((MitData+TString("/Summer13_V1_MC_L2Relative_AK5PF.txt")).Data()); 
    jetCorr->AddCorrectionFromFile((MitData+TString("/Summer13_V1_MC_L3Absolute_AK5PF.txt")).Data()); 
  }
  jetCorr->SetInputName(pubJet->GetOutputName());
  jetCorr->SetCorrectedName("CorrectedJets");    

  mithep::GoodPVFilterMod *goodPvMod = new mithep::GoodPVFilterMod;
  goodPvMod->SetMinVertexNTracks(0);
  goodPvMod->SetMinNDof(4.0);
  goodPvMod->SetMaxAbsZ(24.0);
  goodPvMod->SetMaxRho(2.0);
  goodPvMod->SetIsMC(!isData);
  goodPvMod->SetVertexesName("PrimaryVertexes");

  mithep::MetCorrectionMod *metCorrT0T1Shift = new mithep::MetCorrectionMod;
  metCorrT0T1Shift->SetInputName("PFMet");
  metCorrT0T1Shift->SetJetsName(pubJet->GetOutputName());    
  metCorrT0T1Shift->SetCorrectedJetsName(jetCorr->GetOutputName());    
  metCorrT0T1Shift->SetCorrectedName("PFMetT0T1Shift");   
  metCorrT0T1Shift->ApplyType0(kTRUE);   
  metCorrT0T1Shift->ApplyType1(kTRUE);   
  metCorrT0T1Shift->ApplyShift(kTRUE);   
  metCorrT0T1Shift->IsData(isData);
  metCorrT0T1Shift->SetPrint(kFALSE);

  mithep::MetAnalysisMod* plotMod = new mithep::MetAnalysisMod;
  plotMod->SetMetName(metCorrT0T1Shift->GetOutputName());
  plotMod->SetMetColFromBranch(false);
  plotMod->SetMonoJetCategory(3);

  pubJet->Add(jetCorr);
  jetCorr->Add(goodPvMod);
  goodPvMod->Add(metCorrT0T1Shift);

  ana->SetSuperModule(pubJet);
  ana->AddSuperModule(plotMod);

  ana->Run(false);
}

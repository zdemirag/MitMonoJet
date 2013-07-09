// $Id: runMonoJet.C,v 1.4 2013/07/08 03:45:34 dimatteo Exp $
#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TSystem.h>
#include <TProfile.h>
#include "MitAna/DataUtil/interface/Debug.h"
#include "MitAna/Catalog/interface/Catalog.h"
#include "MitAna/TreeMod/interface/Analysis.h"
#include "MitAna/TreeMod/interface/HLTMod.h"
#include "MitAna/PhysicsMod/interface/RunLumiSelectionMod.h"
#include "MitAna/PhysicsMod/interface/MCProcessSelectionMod.h"
#include "MitAna/PhysicsMod/interface/PublisherMod.h"
#include "MitAna/DataTree/interface/JetCol.h"
#include "MitAna/DataTree/interface/PFJetCol.h"
#include "MitPhysics/Init/interface/ModNames.h"
#include "MitAna/DataTree/interface/Names.h"
#include "MitPhysics/Mods/interface/GoodPVFilterMod.h"
#include "MitPhysics/Mods/interface/MuonIDMod.h"
#include "MitPhysics/Mods/interface/ElectronIDMod.h"
#include "MitPhysics/Mods/interface/ElectronCleaningMod.h"
#include "MitPhysics/Mods/interface/PhotonIDMod.h"
#include "MitPhysics/Mods/interface/PhotonPairSelector.h"
#include "MitPhysics/Mods/interface/PhotonTreeWriter.h"
#include "MitPhysics/Mods/interface/PhotonCleaningMod.h"
#include "MitPhysics/Mods/interface/MergeLeptonsMod.h"
#include "MitPhysics/Mods/interface/JetCorrectionMod.h"
#include "MitPhysics/Mods/interface/PhotonMvaMod.h"
#include "MitPhysics/Mods/interface/MVASystematicsMod.h"
#include "MitPhysics/Mods/interface/SeparatePileUpMod.h"

#include "MitMonoJet/SelMods/interface/MonoJetAnalysisMod.h"
#include "MitMonoJet/Mods/interface/MonoJetTreeWriter.h"

#endif

//--------------------------------------------------------------------------------------------------
void runMonoJet(const char *fileset    = "0000",
				   const char *skim       = "noskim",
				   const char *dataset    = "r12a-pho-j13-v1", 
				   const char *book       = "t2mit/filefi/029",
				   const char *catalogDir = "/home/cmsprod/catalog",
				   const char *outputName = "hgg",
				   int         nEvents    = 1000)
{
  //------------------------------------------------------------------------------------------------
  // some parameters get passed through the environment
  //------------------------------------------------------------------------------------------------
  char json[1024], overlap[1024];
  float overlapCut = -1;
  
  if (gSystem->Getenv("MIT_PROD_JSON"))
    sprintf(json,   "%s",gSystem->Getenv("MIT_PROD_JSON"));
  else {
    sprintf(json, "%s", "~");
    //printf(" JSON file was not properly defined. EXIT!\n");
    //return;
  } 
  
  TString jsonFile = TString("/home/mingyang/cms/json/") + TString(json);
  Bool_t  isData   = ( (jsonFile.CompareTo("/home/mingyang/cms/json/~") != 0) );
  
  if (gSystem->Getenv("MIT_PROD_OVERLAP")) {
    sprintf(overlap,"%s",gSystem->Getenv("MIT_PROD_OVERLAP"));
    if (EOF == sscanf(overlap,"%f",&overlapCut)) {
      printf(" Overlap was not properly defined. EXIT!\n");
      return;
    }
  }
  else {
     sprintf(overlap,"%s", "-1.0");
    //printf(" OVERLAP file was not properly defined. EXIT!\n");
    //return;
  } 

  printf("\n Initialization worked. \n\n");

  //------------------------------------------------------------------------------------------------
  // some global setups
  //------------------------------------------------------------------------------------------------
  using namespace mithep;
  gDebugMask  = Debug::kGeneral;
  gDebugLevel = 3;

  //------------------------------------------------------------------------------------------------
  // set up information
  //------------------------------------------------------------------------------------------------
  RunLumiSelectionMod *runLumiSel = new RunLumiSelectionMod;
  runLumiSel->SetAcceptMC(kTRUE);                          // Monte Carlo events are always accepted

  MCProcessSelectionMod *mcselmod = new MCProcessSelectionMod;
  
  MVASystematicsMod *sysMod = new MVASystematicsMod;
  sysMod->SetMCR9Scale(1.0035, 1.0035);  
  sysMod->SetIsData(isData);
  
  // only select on run- and lumisection numbers when valid json file present
  if ((jsonFile.CompareTo("/home/mingyang/cms/json/~") != 0) &&
      (jsonFile.CompareTo("/home/mingyang/cms/json/-") != 0)   ) {
    runLumiSel->AddJSONFile(jsonFile.Data());
  }
  if ((jsonFile.CompareTo("/home/mingyang/cms/json/-") == 0)   ) {
    printf("\n WARNING -- Looking at data without JSON file: always accept.\n\n");
    runLumiSel->SetAbortIfNotAccepted(kFALSE);   // accept all events if there is no valid JSON file
  }

  printf("\n Run lumi worked. \n\n");

  //------------------------------------------------------------------------------------------------
  // HLT information : [TBF]
  //------------------------------------------------------------------------------------------------
  HLTMod *hltModP = new HLTMod("HLTModP");
  
  //------------------------------------------------------------------------------------------------
  // Run2012
  //------------------------------------------------------------------------------------------------
  hltModP->AddTrigger("HLT_MET120_HBHENoiseCleaned_v*",0,999999); 
  hltModP->AddTrigger("HLT_MonoCentralPFJet80_PFMETnoMu95_NHEF0p95_v*",0,999999); 
  hltModP->AddTrigger("HLT_MonoCentralPFJet80_PFMETnoMu105_NHEF0p95_v*",0,999999);   
    
  hltModP->SetTrigObjsName("MyHltPhotObjs");
  hltModP->SetAbortIfNotAccepted(isData);
  hltModP->SetPrintTable(kFALSE); // set to true to print HLT table
  if(isData == kFALSE){ // ugly, but it works
    hltModP->AddTrigger("HLT_Mu15_v9");
    hltModP->AddTrigger("!HLT_Mu15_v9");
    hltModP->AddTrigger("HLT_Mu15_v2");
    hltModP->AddTrigger("!HLT_Mu15_v2");
    hltModP->AddTrigger("HLT_Mu9");
    hltModP->AddTrigger("!HLT_Mu9");
    hltModP->AddTrigger("HLT_Mu12_v13");
    hltModP->AddTrigger("!HLT_Mu12_v13");
    hltModP->AddTrigger("HLT_Mu12_v14");
    hltModP->AddTrigger("!HLT_Mu12_v14");
    hltModP->AddTrigger("HLT_Mu12_v16");
    hltModP->AddTrigger("!HLT_Mu12_v16");
    // // MC signal triggers
    // hltModP->AddTrigger("HLT_MET120_HBHENoiseCleaned_v*");
    // hltModP->AddTrigger("HLT_MonoCentralPFJet80_PFMETnoMu95_NHEF0p95_v*");
  }

  //------------------------------------------------------------------------------------------------
  // split pfcandidates to PFPU and PFnoPU
  //------------------------------------------------------------------------------------------------
  SeparatePileUpMod* SepPUMod = new SeparatePileUpMod;
  //  SepPUMod->SetUseAllVerteces(kFALSE);
  // SepPUMod->SetVertexName("OutVtxCiC");
  SepPUMod->SetPFNoPileUpName("pfnopileupcands");
  SepPUMod->SetPFPileUpName("pfpileupcands");
  SepPUMod->SetCheckClosestZVertex(kFALSE);
  
  //------------------------------------------------------------------------------------------------
  // select events with a good primary vertex
  //------------------------------------------------------------------------------------------------
  GoodPVFilterMod *goodPVFilterMod = new GoodPVFilterMod;
  goodPVFilterMod->SetMinVertexNTracks(0);
  goodPVFilterMod->SetMinNDof         (4.0);
  goodPVFilterMod->SetMaxAbsZ         (24.0);
  goodPVFilterMod->SetMaxRho          (2.0);
  goodPVFilterMod->SetAbortIfNotAccepted(kFALSE);
  goodPVFilterMod->SetIsMC(!isData);
  
  PublisherMod<PFJet,Jet> *pubJet = new PublisherMod<PFJet,Jet>("JetPub");
  pubJet->SetInputName("AKt5PFJets");
  pubJet->SetOutputName("PubAKt5PFJets");
  
  PublisherMod<PFJet,Jet> *pubJetOpen = new PublisherMod<PFJet,Jet>("JetPubOpen");
  pubJetOpen->SetInputName("AKt5PFJets");
  pubJetOpen->SetOutputName("PubAKt5PFJetsOpen");  
  
  JetCorrectionMod *jetCorr = new JetCorrectionMod;
  if(isData){ 
    jetCorr->AddCorrectionFromFile(std::string((gSystem->Getenv("CMSSW_BASE") + TString("/src/MitPhysics/data/Summer12_V7_DATA_L1FastJet_AK5PF.txt")).Data())); 
    jetCorr->AddCorrectionFromFile(std::string((gSystem->Getenv("CMSSW_BASE") + TString("/src/MitPhysics/data/Summer12_V7_DATA_L2Relative_AK5PF.txt")).Data())); 
    jetCorr->AddCorrectionFromFile(std::string((gSystem->Getenv("CMSSW_BASE") + TString("/src/MitPhysics/data/Summer12_V7_DATA_L3Absolute_AK5PF.txt")).Data())); 
    jetCorr->AddCorrectionFromFile(std::string((gSystem->Getenv("CMSSW_BASE") + TString("/src/MitPhysics/data/Summer12_V7_DATA_L2L3Residual_AK5PF.txt")).Data())); 
  }
  else {
    jetCorr->AddCorrectionFromFile(std::string((gSystem->Getenv("CMSSW_BASE") + TString("/src/MitPhysics/data/Summer12_V7_MC_L1FastJet_AK5PF.txt")).Data())); 
    jetCorr->AddCorrectionFromFile(std::string((gSystem->Getenv("CMSSW_BASE") + TString("/src/MitPhysics/data/Summer12_V7_MC_L2Relative_AK5PF.txt")).Data())); 
    jetCorr->AddCorrectionFromFile(std::string((gSystem->Getenv("CMSSW_BASE") + TString("/src/MitPhysics/data/Summer12_V7_MC_L3Absolute_AK5PF.txt")).Data())); 
  }
  jetCorr->SetInputName(pubJet->GetOutputName());
  jetCorr->SetCorrectedName("CorrectedJets");    
  
  
  Bool_t excludedoubleprompt = kFALSE;
  if (TString(dataset).Contains("-pj")) {
    mcselmod->ExcludeProcess(18);
    mcselmod->ExcludeProcess(114);
    excludedoubleprompt = kTRUE;
  }
  
  if (TString(dataset).Contains("-qcd2em") || TString(dataset).Contains("-qcd-2em")) {
    excludedoubleprompt = kTRUE;
  }
   
  //------------------------------------------------------------------------------------------------
  // select events with jet+MET
  //------------------------------------------------------------------------------------------------
  MonoJetAnalysisMod         *jetplusmet = new MonoJetAnalysisMod("MonoJetSelector");
  jetplusmet->SetJetsName(jetCorr->GetOutputName());
  jetplusmet->SetJetsFromBranch(kFALSE);
  jetplusmet->SetMinNumJets(1);
  jetplusmet->SetMinJetEt(40);
  jetplusmet->SetMaxJetEta(2.4);
  jetplusmet->SetMinMetEt(40);

  MonoJetTreeWriter *jetplusmettree = new MonoJetTreeWriter("MonoJetTreeWriter");
  jetplusmettree->SetPhotonsFromBranch(kTRUE);
  jetplusmettree->SetEnableJets(kTRUE);
  jetplusmettree->SetPFJetsFromBranch(kFALSE);
  jetplusmettree->SetPFJetName(jetCorr->GetOutputName());
  jetplusmettree->SetIsData(isData);
 
  //------------------------------------------------------------------------------------------------
  // making analysis chain
  //------------------------------------------------------------------------------------------------
  // this is how it always starts
  runLumiSel      ->Add(mcselmod);
    
  if (TString(dataset).Contains("-h")) {    
    mcselmod        ->Add(sysMod);
  }
    
  if (!TString(dataset).Contains("meridiani")) {      
 //   mcselmod         ->Add(hltModP);
 //   hltModP         ->Add(goodPVFilterMod);
    mcselmod->Add(goodPVFilterMod);
  }
  else {
    mcselmod->Add(goodPVFilterMod);
  }
  //hltModE         ->Add(goodPVFilterModE);
  //hltModES        ->Add(goodPVFilterModES);
  
  
  // simple object id modules
  goodPVFilterMod  ->Add(pubJet);
  pubJet           ->Add(jetCorr);
  jetCorr          ->Add(SepPUMod); 
   
  // Test jet+met ana
  SepPUMod         ->Add(jetplusmet);
  jetplusmet       ->Add(jetplusmettree);

  //TFile::SetOpenTimeout(0);
  //TFile::SetCacheFileDir("./rootfilecache",kTRUE,kTRUE);
  TFile::SetReadaheadSize(128*1024*1024);
  
  //------------------------------------------------------------------------------------------------
  // setup analysis
  //------------------------------------------------------------------------------------------------
  Analysis *ana = new Analysis;
  ana->SetUseHLT(kTRUE);
  ana->SetKeepHierarchy(kTRUE);
  ana->SetSuperModule(runLumiSel);
  ana->SetPrintScale(100);
  if (nEvents >= 0)
    ana->SetProcessNEvents(nEvents);
  
  //------------------------------------------------------------------------------------------------
  // organize input
  //------------------------------------------------------------------------------------------------
  Catalog *c = new Catalog(catalogDir);
  TString skimdataset = TString(dataset)+TString("/") +TString(skim);
  Dataset *d = NULL;
  TString bookstr = book;
  //if (TString(dataset).Contains("s11-h")) bookstr.ReplaceAll("local","t2mit");
  if (TString(skim).CompareTo("noskim") == 0)
    d = c->FindDataset(bookstr,dataset,fileset,true);
  else 
    d = c->FindDataset(bookstr,skimdataset.Data(),fileset,true);
  ana->AddDataset(d);
  // // only for testing
  //ana->AddFile("/mnt/hadoop/cmsprod/filefi/029/s12-h125gg-gf-v7a/D0515547-68FA-E111-B849-002618FDA262.root");

  //------------------------------------------------------------------------------------------------
  // organize output
  //------------------------------------------------------------------------------------------------
  TString rootFile = TString(outputName);
  rootFile += TString("_") + TString(dataset) + TString("_") + TString(skim);
  if (TString(fileset) != TString(""))
    rootFile += TString("_") + TString(fileset);
  rootFile += TString(".root");
  ana->SetOutputName(rootFile.Data());
  //ana->SetCacheSize(64*1024*1024);
  ana->SetCacheSize(0);
  
  //------------------------------------------------------------------------------------------------
  // Say what we are doing
  //------------------------------------------------------------------------------------------------
  printf("\n==== PARAMETER SUMMARY FOR THIS JOB ====\n");
  printf("\n JSON file: %s\n  and overlap cut: %f (%s)\n",jsonFile.Data(),overlapCut,overlap);
  printf("\n Rely on Catalog: %s\n",catalogDir);
  printf("  -> Book: %s  Dataset: %s  Skim: %s  Fileset: %s <-\n",book,dataset,skim,fileset);
  printf("\n Root output: %s\n\n",rootFile.Data());  
  printf("\n========================================\n");

  //------------------------------------------------------------------------------------------------
  // run the analysis after successful initialisation
  //------------------------------------------------------------------------------------------------
  ana->Run(!gROOT->IsBatch());

  return;
}

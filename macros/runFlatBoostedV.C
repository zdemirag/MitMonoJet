#include <TSystem.h>
#include <TProfile.h>
#include "MitCommon/Utils/interface/Utils.h"
#include "MitAna/DataUtil/interface/Debug.h"
#include "MitAna/Catalog/interface/Catalog.h"
#include "MitAna/TreeMod/interface/Analysis.h"
#include "MitAna/TreeMod/interface/HLTMod.h"
#include "MitAna/PhysicsMod/interface/RunLumiSelectionMod.h"
#include "MitAna/PhysicsMod/interface/MCProcessSelectionMod.h"
#include "MitAna/PhysicsMod/interface/PublisherMod.h"
#include "MitAna/DataTree/interface/Names.h"
#include "MitPhysics/Init/interface/ModNames.h"
#include "MitPhysics/Utils/interface/RhoUtilities.h"
#include "MitPhysics/Mods/interface/GeneratorMod.h"
#include "MitPhysics/Mods/interface/GoodPVFilterMod.h"
#include "MitMonoJet/Mods/interface/DMSTreeWriter.h"

TString getCatalogDir(const char* dir);
TString getJsonFile(const char* dir);

//--------------------------------------------------------------------------------------------------
void runFlatBoostedV(const char *fileset    = "0000",
                     const char *skim       = "noskim",
                     const char *dataset    = "s12-dmmjet-avd_m1-v7a",     
                     const char *book       = "t2mit/filefi/032",
                     const char *catalogDir = "/home/cmsprod/catalog",
                     const char *outputName = "boostedv",
                     int         nEvents    = 100)
{
  //------------------------------------------------------------------------------------------------
  // some parameters get passed through the environment
  //------------------------------------------------------------------------------------------------
  TString cataDir  = getCatalogDir(catalogDir);
  TString mitData  = Utils::GetEnv("MIT_DATA");
  TString json     = Utils::GetEnv("MIT_PROD_JSON");
  TString jsonFile = getJsonFile("/home/cmsprod/cms/json");
  Bool_t  isData   = (json.CompareTo("~") != 0);
  printf("\n Initialization worked. Data?: %d\n\n",isData);

  //------------------------------------------------------------------------------------------------
  // some global setups
  //------------------------------------------------------------------------------------------------
  // debugging config
  using namespace mithep;
  gDebugMask  = (Debug::EDebugMask) (Debug::kGeneral | Debug::kTreeIO);
  gDebugLevel = 3;


  // Caching and how
  Int_t local = 1, cacher = 1;

  // local =   0 - as is,
  //           1 - /mt/hadoop  (MIT:SmartCache - preload one-by-one)
  //           2 - /mnt/hadoop (MIT:SmartCache - preload complete fileset)
  //           3 - ./          (xrdcp          - preload one-by-one)
  // cacher =  0 - no file by file caching
  //           1 - file by file caching on

  //------------------------------------------------------------------------------------------------
  // set up information for master module
  //------------------------------------------------------------------------------------------------
  RunLumiSelectionMod *runLumiSel = new RunLumiSelectionMod;
  runLumiSel->SetAcceptMC(kTRUE);                          // Monte Carlo events are always accepted
  
  // only select on run- and lumisection numbers when valid json file present
  if (json.CompareTo("~") != 0 && json.CompareTo("-") != 0) {
    printf(" runBoostedV() - adding jsonFile: %s\n",jsonFile.Data());
    runLumiSel->AddJSONFile(jsonFile.Data());
  }
  if (json.CompareTo("-") == 0) {
    printf("\n WARNING -- Looking at data without JSON file: always accept.\n\n");
    runLumiSel->SetAbortIfNotAccepted(kFALSE);   // accept all events if there is no valid JSON file
  }
  printf("\n Run lumi worked. \n\n");

  //------------------------------------------------------------------------------------------------
  // setup analysis
  //------------------------------------------------------------------------------------------------
  Analysis *ana = new Analysis;
  ana->SetUseCacher(cacher);
  ana->SetUseHLT(kFALSE);
  ana->SetKeepHierarchy(kTRUE);
  ana->SetSuperModule(runLumiSel);
  ana->SetPrintScale(100);
  if (nEvents >= 0)
    ana->SetProcessNEvents(nEvents);

  //------------------------------------------------------------------------------------------------
  // organize input
  //------------------------------------------------------------------------------------------------
  Catalog *c = new Catalog(cataDir.Data());
  TString skimdataset = TString(dataset)+TString("/") +TString(skim);
  Dataset *d = NULL;
  TString bookstr = book;
  //if (TString(skim).CompareTo("noskim") == 0)
  //  d = c->FindDataset(bookstr,dataset,fileset,local);
  //else 
  //  d = c->FindDataset(bookstr,skimdataset.Data(),fileset,local);
  ana->AddFile("/scratch1/dimatteo/cmssw/033/CMSSW_5_3_15/src/MitMonoJet/macros/boostedv_s12-dmmjet-avd_m1-v7a_noskim_0000_ntuple_000.root");

  //------------------------------------------------------------------------------------------------
  // organize output
  //------------------------------------------------------------------------------------------------
  TString rootFile = TString(outputName);
  rootFile += TString("_") + TString(dataset) + TString("_") + TString(skim);
  if (TString(fileset) != TString(""))
    rootFile += TString("_") + TString(fileset);
  TString ntupleFile = rootFile + TString("_ntuple");
  rootFile += TString(".root");
  ana->SetOutputName(rootFile.Data());
  ana->SetCacheSize(0);

  //------------------------------------------------------------------------------------------------
  // HLT information
  //------------------------------------------------------------------------------------------------
  //HLTMod *hltModP = new HLTMod("HLTModP");
  //hltModP->SetBitsName("HLTBits");
  //hltModP->SetTrigObjsName("HltObjsMonoJet");
  //hltModP->SetAbortIfNotAccepted(isData);
  //hltModP->SetPrintTable(kFALSE);

  //// monojet triggers
  //const int nMjtTrigs = 12;
  //TString monoJetTriggers[nMjtTrigs] = { "HLT_MonoCentralPFJet80_PFMETnoMu105_NHEF0p95_v4",
                                         //"HLT_MonoCentralPFJet80_PFMETnoMu105_NHEF0p95_v3",
                                         //"HLT_MonoCentralPFJet80_PFMETnoMu105_NHEF0p95_v1",
                                         //"HLT_MonoCentralPFJet80_PFMETnoMu95_NHEF0p95_v5",
                                         //"HLT_MonoCentralPFJet80_PFMETnoMu95_NHEF0p95_v4",
                                         //"HLT_MonoCentralPFJet80_PFMETnoMu95_NHEF0p95_v3",
                                         //"HLT_MonoCentralPFJet80_PFMETnoMu95_NHEF0p95_v2",
                                         //"HLT_MET120_HBHENoiseCleaned_v6",
                                         //"HLT_MET120_HBHENoiseCleaned_v5",
                                         //"HLT_MET120_HBHENoiseCleaned_v4",
                                         //"HLT_MET120_HBHENoiseCleaned_v3",
                                         //"HLT_MET120_HBHENoiseCleaned_v2" };

  //for (int i=0; i<nMjtTrigs; i++)
    //hltModP->AddTrigger(TString("!+"+monoJetTriggers[i]),0,999999);

  //------------------------------------------------------------------------------------------------
  // prepare the tree writer 
  //------------------------------------------------------------------------------------------------
  DMSTreeWriter *treeWriter = new DMSTreeWriter();
  treeWriter->SetIsData(kFALSE);
  treeWriter->SetMetName("SkmPFMetT0T1Shift");
  treeWriter->SetMetMVAName("PFMetMVA");
  treeWriter->SetPhotonsName("SkmCleanPhotons");
  treeWriter->SetElectronsName("SkmCleanElectrons");
  treeWriter->SetMuonsName("SkmHWWMuons");
  treeWriter->SetTausName("SkmCleanPFTaus");
  treeWriter->SetJetsName("SkmCleanJets");
  treeWriter->SetTriggerObjectsName("SkmHltObjsMonoJet");
  
  //------------------------------------------------------------------------------------------------
  // making analysis chain
  //------------------------------------------------------------------------------------------------
  runLumiSel               ->Add(treeWriter);
  
  //------------------------------------------------------------------------------------------------
  // Say what we are doing
  //------------------------------------------------------------------------------------------------
  printf("\n==== PARAMETER SUMMARY FOR THIS JOB ====\n");
  printf("\n JSON file: %s\n",jsonFile.Data());
  printf("\n Rely on Catalog: %s\n",cataDir.Data());
  printf("  -> Book: %s  Dataset: %s  Skim: %s  Fileset: %s <-\n",book,dataset,skim,fileset);
  printf("\n Root output:   %s\n",rootFile.Data());  
  printf("\n Ntuple output: %s\n\n",(ntupleFile + TString(".root")).Data());  
  printf("\n========================================\n");

  //------------------------------------------------------------------------------------------------
  // run the analysis after successful initialisation
  //------------------------------------------------------------------------------------------------
  ana->Run(!gROOT->IsBatch());

  return;
}

//--------------------------------------------------------------------------------------------------
TString getCatalogDir(const char* dir)
{
  TString cataDir = TString("./catalog");
  Long_t *id=0,*size=0,*flags=0,*mt=0;

  printf(" Try local catalog first: %s\n",cataDir.Data());
  if (gSystem->GetPathInfo(cataDir.Data(),id,size,flags,mt) != 0) {
    cataDir = TString(dir);
    if (gSystem->GetPathInfo(cataDir.Data(),id,size,flags,mt) != 0) {
      printf(" Requested local (./catalog) and specified catalog do not exist. EXIT!\n");
      return TString("");
    }
  }
  else {
    printf(" Local catalog exists: %s using this one.\n",cataDir.Data()); 
  }

  return cataDir;
}

//--------------------------------------------------------------------------------------------------
TString getJsonFile(const char* dir)
{
  TString jsonDir  = TString("./json");
  TString json     = Utils::GetEnv("MIT_PROD_JSON");
  Long_t *id=0,*size=0,*flags=0,*mt=0;

  printf(" Try local json first: %s\n",jsonDir.Data());
  if (gSystem->GetPathInfo(jsonDir.Data(),id,size,flags,mt) != 0) {
    jsonDir = TString(dir);
    if (gSystem->GetPathInfo(jsonDir.Data(),id,size,flags,mt) != 0) {
      printf(" Requested local (./json) and specified json directory do not exist. EXIT!\n");
      return TString("");
    }
  }
  else {
    printf(" Local json directory exists: %s using this one.\n",jsonDir.Data()); 
  }

  // Construct the full file name
  TString jsonFile = jsonDir + TString("/") + json;
  if (gSystem->GetPathInfo(jsonFile.Data(),id,size,flags,mt) != 0) {
    printf(" Requested jsonfile (%s) does not exist. EXIT!\n",jsonFile.Data());
    return TString("");
  }
  else {
    printf(" Requested jsonfile (%s) exist. Moving on now!\n",jsonFile.Data());
  }

  return jsonFile;
}

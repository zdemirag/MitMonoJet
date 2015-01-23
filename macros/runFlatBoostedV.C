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
                     const char *dataset    = "s12-ttj-v1-v7a",     
                     const char *book       = "t2mit/filefi/032",
                     const char *catalogDir = "/home/cmsprod/catalog",
                     const char *outputName = "boostedv-v9",
                     int         nEvents    = -1)
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
  Int_t cacher = 0;

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
  TString skimdataset = TString(dataset)+TString("/") +TString(skim);
  TString bookstr = book;
  TString inputFileList;
  inputFileList = "inputBavanti.txt";
  ana->AddFiles(inputFileList,-1);

  TString inputPUFile = Utils::GetEnv("MIT_PROD_HIST") + "/" + TString(outputName) + "/merged/";
  inputPUFile += TString(outputName) + TString("_") +  TString(dataset);
  inputPUFile += "_noskim.root";  

  //------------------------------------------------------------------------------------------------
  // organize output
  //------------------------------------------------------------------------------------------------
  TString rootFile = TString(outputName);
  rootFile += TString("_") + TString(dataset) + TString("_") + TString(skim);
  rootFile += TString("_flatntuple") + TString(".root");
  ana->SetOutputName(rootFile.Data());
  ana->SetCacheSize(0);

  //------------------------------------------------------------------------------------------------
  // prepare the tree writer 
  //------------------------------------------------------------------------------------------------
  DMSTreeWriter *treeWriter = new DMSTreeWriter();
  treeWriter->SetIsData(isData);
  treeWriter->SetPhotonsName("XsPhotons");
  treeWriter->SetElectronsName("XsElectrons");
  treeWriter->SetMuonsName("XsMuons");
  treeWriter->SetTausName("XsTaus");
  treeWriter->SetJetsName("XlJets");
  treeWriter->SetTriggerObjectsName("SkmHltObjsMonoJet");
  treeWriter->SetInPUHistoFileName(inputPUFile);
  treeWriter->SetTargetPUHistoFileName("/home/dimatteo/cms/external/MyDataPileupHistogram.root");
  
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
  printf("\n========================================\n");

  //------------------------------------------------------------------------------------------------
  // run the analysis after successful initialisation
  //------------------------------------------------------------------------------------------------
  ana->Run(!gROOT->IsBatch());

  //------------------------------------------------------------------------------------------------
  // copy the file to its final destination
  //------------------------------------------------------------------------------------------------
  //gSystem->Exec("mv " + rootFile + TString(" /scratch4/dimatteo/cms/hist/") + TString(outputName) + TString("/merged/")); 

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

#if !defined(__CINT__) || defined(__MAKECINT__)
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
#include "MitAna/DataTree/interface/JetCol.h"
#include "MitAna/DataTree/interface/PFJetCol.h"
#include "MitPhysics/Init/interface/ModNames.h"
#include "MitAna/DataTree/interface/Names.h"
#include "MitPhysics/Mods/interface/GeneratorMod.h"
#include "MitPhysics/Mods/interface/GoodPVFilterMod.h"
#include "MitPhysics/Mods/interface/MuonIDMod.h"
#include "MitPhysics/Mods/interface/ElectronIDMod.h"
#include "MitPhysics/Mods/interface/ElectronCleaningMod.h"
#include "MitPhysics/Mods/interface/PFTauIDMod.h"
#include "MitPhysics/Mods/interface/PFTauCleaningMod.h"
#include "MitPhysics/Mods/interface/PhotonIDMod.h"
#include "MitPhysics/Mods/interface/PhotonTreeWriter.h"
#include "MitPhysics/Mods/interface/PhotonCleaningMod.h"
#include "MitPhysics/Mods/interface/MergeLeptonsMod.h"
#include "MitPhysics/Mods/interface/JetCorrectionMod.h"
#include "MitPhysics/Mods/interface/MetCorrectionMod.h"
#include "MitPhysics/Mods/interface/PhotonMvaMod.h"
#include "MitPhysics/Mods/interface/MVASystematicsMod.h"
#include "MitPhysics/Mods/interface/SeparatePileUpMod.h"
#include "MitPhysics/Mods/interface/JetIDMod.h"
#include "MitPhysics/Mods/interface/JetCleaningMod.h"
#include "MitMonoJet/SelMods/interface/MonoJetAnalysisMod.h"
#include "MitMonoJet/Mods/interface/MonoJetTreeWriter.h"
#endif

TString getCatalogDir(const char* dir);
TString getJsonFile(const char* dir);

//--------------------------------------------------------------------------------------------------
void runMonoJet(const char *fileset    = "0000",
                const char *skim       = "noskim",
                const char *dataset    = "s12-wjets-1-v7a",
                const char *book       = "t2mit/filefi/032",
                const char *catalogDir = "/home/cmsprod/catalog",
                const char *outputName = "MonoJet_August13",
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

  std::cout<<"*********** Is data?? **********"<<isData<<std::endl;

  //------------------------------------------------------------------------------------------------
  // some global setups
  //------------------------------------------------------------------------------------------------
  using namespace mithep;
  gDebugMask  = (Debug::EDebugMask) (Debug::kGeneral | Debug::kTreeIO);
  gDebugLevel = 3;

  // Caching and how
  Int_t local = 1, cacher = 1;

  // local =   0 - as is,
  //           1 - /mnt/hadoop (MIT:SmartCache - preload one-by-one)
  //           2 - /mnt/hadoop (MIT:SmartCache - preload complete fileset)
  //           3 - ./          (xrdcp          - preload one-by-one)
  // cacher =  0 - no file by file caching
  //           1 - file by file caching on

  //------------------------------------------------------------------------------------------------
  // set up information
  //------------------------------------------------------------------------------------------------
  RunLumiSelectionMod *runLumiSel = new RunLumiSelectionMod;
  runLumiSel->SetAcceptMC(!isData);
  runLumiSel->SetAbortIfNotAccepted(kFALSE);   // accept all events if there is no valid JSON file

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

  // Generator info
  GeneratorMod *generatorMod = new GeneratorMod;
  generatorMod->SetPrintDebug(kFALSE);
  generatorMod->SetPtLeptonMin(0.0);
  generatorMod->SetEtaLeptonMax(2.7);
  generatorMod->SetPtPhotonMin(0.0);
  generatorMod->SetEtaPhotonMax(2.7);
  generatorMod->SetPtRadPhotonMin(0.0);
  generatorMod->SetEtaRadPhotonMax(2.7);
  generatorMod->SetIsData(isData);
  generatorMod->SetFillHist(!isData);
  generatorMod->SetApplyISRFilter(kFALSE);
  generatorMod->SetApplyVVFilter(kFALSE);
  generatorMod->SetApplyVGFilter(kFALSE);
  generatorMod->SetFilterBTEvents(kFALSE);

  //-----------------------------------------------------------------------------------------------------------
  // HLT information : trigger not applied (neither for data nor for MC, store info to apply selection offline
  //-----------------------------------------------------------------------------------------------------------
  HLTMod *hltModP = new HLTMod("HLTModP");

  // monojet triggers
  const int nMjtTrigs = 12;
  TString monoJetTriggers[nMjtTrigs] = { "HLT_MonoCentralPFJet80_PFMETnoMu105_NHEF0p95_v4",
                                         "HLT_MonoCentralPFJet80_PFMETnoMu105_NHEF0p95_v3",
                                         "HLT_MonoCentralPFJet80_PFMETnoMu105_NHEF0p95_v1",
                                         "HLT_MonoCentralPFJet80_PFMETnoMu95_NHEF0p95_v5",
                                         "HLT_MonoCentralPFJet80_PFMETnoMu95_NHEF0p95_v4",
                                         "HLT_MonoCentralPFJet80_PFMETnoMu95_NHEF0p95_v3",
                                         "HLT_MonoCentralPFJet80_PFMETnoMu95_NHEF0p95_v2",
                                         "HLT_MET120_HBHENoiseCleaned_v6",
                                         "HLT_MET120_HBHENoiseCleaned_v5",
                                         "HLT_MET120_HBHENoiseCleaned_v4",
                                         "HLT_MET120_HBHENoiseCleaned_v3",
                                         "HLT_MET120_HBHENoiseCleaned_v2" };

  for (int i=0; i<nMjtTrigs; i++)
    hltModP->AddTrigger(TString("!+"+monoJetTriggers[i]),0,999999);

  // VBF triggers
  const int nVbfTrigs = 7;
  TString vbfTriggers[nVbfTrigs] = { "HLT_DiPFJet40_PFMETnoMu65_MJJ800VBF_AllJets_v9",
                                     "HLT_DiPFJet40_PFMETnoMu65_MJJ800VBF_AllJets_v8",
                                     "HLT_DiPFJet40_PFMETnoMu65_MJJ800VBF_AllJets_v6",
                                     "HLT_DiPFJet40_PFMETnoMu65_MJJ800VBF_AllJets_v5",
                                     "HLT_DiPFJet40_PFMETnoMu65_MJJ800VBF_AllJets_v4",
                                     "HLT_DiPFJet40_PFMETnoMu65_MJJ800VBF_AllJets_v3",
                                     "HLT_DiPFJet40_PFMETnoMu65_MJJ800VBF_AllJets_v2" };

  for (int i=0; i<nVbfTrigs; i++)
    hltModP->AddTrigger((TString("!+")+vbfTriggers[i]).Data(),0,999999);

  hltModP->SetBitsName("HLTBits");
  hltModP->SetTrigObjsName("MyHltPhotObjs");
  hltModP->SetAbortIfNotAccepted(isData);
  hltModP->SetPrintTable(kFALSE);

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
  goodPVFilterMod->SetIsMC(!isData);
  goodPVFilterMod->SetVertexesName("PrimaryVertexes");

  //------------------------------------------------------------------------------------------------
  // object id and cleaning sequence
  //------------------------------------------------------------------------------------------------
  //-----------------------------------
  // Lepton Selection
  //-----------------------------------
  ElectronIDMod* eleIdMod = new ElectronIDMod;
  eleIdMod->SetPtMin(10.);
  eleIdMod->SetEtaMax(2.5);
  eleIdMod->SetApplyEcalFiducial(kTRUE);
  eleIdMod->SetIDType("CustomLoose");
  eleIdMod->SetIsoType("PFIso_HggLeptonTag2012HCP");
  eleIdMod->SetPFNoPileUpName("pfnopileupcands");
  eleIdMod->SetApplyConversionFilterType1(kTRUE);
  eleIdMod->SetApplyConversionFilterType2(kFALSE);
  eleIdMod->SetChargeFilter(kFALSE);
  eleIdMod->SetApplyD0Cut(kTRUE);
  eleIdMod->SetApplyDZCut(kTRUE);
  eleIdMod->SetWhichVertex(0);
  eleIdMod->SetNExpectedHitsInnerCut(0);
  eleIdMod->SetGoodElectronsName("GoodElectronsBS");
  eleIdMod->SetRhoType(RhoUtilities::CMS_RHO_RHOKT6PFJETS);

  MuonIDMod *muonIdWW = new MuonIDMod;
  muonIdWW->SetOutputName("HWWMuons");
  muonIdWW->SetIntRadius(0.0);
  muonIdWW->SetClassType("GlobalTracker");
  muonIdWW->SetIDType("WWMuIdV4");
  muonIdWW->SetIsoType("IsoRingsV0_BDTG_Iso");
  muonIdWW->SetApplyD0Cut(kTRUE);
  muonIdWW->SetApplyDZCut(kTRUE);
  muonIdWW->SetWhichVertex(0);
  muonIdWW->SetRhoType(RhoUtilities::CMS_RHO_RHOKT6PFJETS);
  muonIdWW->SetPtMin(10.);
  muonIdWW->SetEtaCut(2.4);

  MuonIDMod* muonIdIso = new MuonIDMod;
  muonIdIso->SetOutputName("IsolatedPOGMuons");
  muonIdIso->SetClassType("GlobalorTracker");
  muonIdIso->SetIDType("NoId");
  muonIdIso->SetApplyD0Cut(true);
  muonIdIso->SetD0Cut(0.2);
  muonIdIso->SetApplyDZCut(true);
  muonIdIso->SetDZCut(0.5);
  muonIdIso->SetIsoType("PFIsoBetaPUCorrected"); //h
  muonIdIso->SetPFNoPileUpName("pfnopileupcands");
  muonIdIso->SetPFPileUpName("pfpileupcands");
  muonIdIso->SetPtMin(10.);
  muonIdIso->SetEtaCut(2.4);

  MuonIDMod *muonIdPOG = new MuonIDMod;
  muonIdPOG->SetOutputName("POGMuons");
  muonIdPOG->SetClassType("GlobalTracker");
  muonIdPOG->SetIDType("NoId");
  muonIdPOG->SetApplyD0Cut(true);
  muonIdPOG->SetD0Cut(0.2);
  muonIdPOG->SetApplyDZCut(true);
  muonIdPOG->SetDZCut(0.5);
  muonIdPOG->SetIsoType("NoIso");
  muonIdPOG->SetPtMin(10.);
  muonIdPOG->SetEtaCut(2.4);


  MuonIDMod *muonIdLoose = new MuonIDMod;
  muonIdLoose->SetOutputName("LooseMuons");
  muonIdLoose->SetClassType("All");
  muonIdLoose->SetIDType("NoId");
  muonIdLoose->SetIsoType("NoIso");
  muonIdLoose->SetApplyD0Cut(false);
  muonIdLoose->SetApplyDZCut(false);
  muonIdLoose->SetPtMin(10.);
  muonIdLoose->SetEtaCut(2.4);

  //MuonIDMod *muonId = muonIdPOG;
  MuonIDMod *muonId = muonIdLoose;

  ElectronCleaningMod *electronCleaning = new ElectronCleaningMod;
  electronCleaning->SetCleanMuonsName(muonId->GetOutputName());
  electronCleaning->SetGoodElectronsName(eleIdMod->GetOutputName());
  electronCleaning->SetCleanElectronsName("CleanElectrons");

  MergeLeptonsMod *merger = new MergeLeptonsMod;
  merger->SetMuonsName(muonId->GetOutputName());
  merger->SetElectronsName(electronCleaning->GetOutputName());
  merger->SetMergedName("MergedLeptons");

//-----------------------------------
  // Photon Regression + ID
  //-----------------------------------
  PhotonMvaMod *photreg = new PhotonMvaMod;
  photreg->SetRegressionVersion(3);
  photreg->SetRegressionWeights(std::string(
    (gSystem->Getenv("MIT_DATA") + TString("/gbrv3ph_52x.root")).Data()
    ));
  photreg->SetOutputName("GoodPhotonsRegr");
  photreg->SetApplyShowerRescaling(kTRUE);
  photreg->SetMinNumPhotons(0);
  photreg->SetIsData(isData);

  PhotonIDMod *photonIDMod = new PhotonIDMod;
  photonIDMod->SetPtMin(15.0);
  photonIDMod->SetOutputName("GoodPhotons");
  photonIDMod->SetIDType("MITMVAId");
  photonIDMod->SetBdtCutBarrel(0.02);
  photonIDMod->SetBdtCutEndcap(0.1);
  photonIDMod->SetIdMVAType("2013FinalIdMVA_8TeV");
  photonIDMod->SetApplyElectronVeto(kTRUE);
  photonIDMod->SetApplyFiduciality(kTRUE);
  photonIDMod->SetIsData(isData);
  photonIDMod->SetPhotonsFromBranch(kFALSE);
  photonIDMod->SetInputName(photreg->GetOutputName());

  PFTauIDMod *pftauIDMod = new PFTauIDMod;
  pftauIDMod->SetPFTausName("HPSTaus");
  pftauIDMod->SetIsLooseId(kFALSE);
  pftauIDMod->SetIsHPSSel(kTRUE); // to get >= 5_3_14 samples running

  PhotonCleaningMod *photonCleaningMod = new PhotonCleaningMod;
  photonCleaningMod->SetCleanElectronsName(electronCleaning->GetOutputName());
  photonCleaningMod->SetGoodPhotonsName(photonIDMod->GetOutputName());
  photonCleaningMod->SetCleanPhotonsName("CleanPhotons");

  PFTauCleaningMod *pftauCleaningMod = new PFTauCleaningMod;
  pftauCleaningMod->SetGoodPFTausName(pftauIDMod->GetGoodPFTausName());
  pftauCleaningMod->SetCleanMuonsName(muonId->GetOutputName());

  PublisherMod<PFJet,Jet> *pubJet = new PublisherMod<PFJet,Jet>("JetPub");
  pubJet->SetInputName("AKt5PFJets");
  pubJet->SetOutputName("PubAKt5PFJets");

  JetCorrectionMod *jetCorr = new JetCorrectionMod;
  if(isData){
    jetCorr->AddCorrectionFromFile(std::string((gSystem->Getenv("CMSSW_BASE") + TString("/src/MitPhysics/data/Summer13_V1_DATA_L1FastJet_AK5PF.txt")).Data()));
    jetCorr->AddCorrectionFromFile(std::string((gSystem->Getenv("CMSSW_BASE") + TString("/src/MitPhysics/data/Summer13_V1_DATA_L2Relative_AK5PF.txt")).Data()));
    jetCorr->AddCorrectionFromFile(std::string((gSystem->Getenv("CMSSW_BASE") + TString("/src/MitPhysics/data/Summer13_V1_DATA_L3Absolute_AK5PF.txt")).Data()));
    jetCorr->AddCorrectionFromFile(std::string((gSystem->Getenv("CMSSW_BASE") + TString("/src/MitPhysics/data/Summer13_V1_DATA_L2L3Residual_AK5PF.txt")).Data()));
  }
  else {
    jetCorr->AddCorrectionFromFile(std::string((gSystem->Getenv("CMSSW_BASE") + TString("/src/MitPhysics/data/Summer13_V1_MC_L1FastJet_AK5PF.txt")).Data()));
    jetCorr->AddCorrectionFromFile(std::string((gSystem->Getenv("CMSSW_BASE") + TString("/src/MitPhysics/data/Summer13_V1_MC_L2Relative_AK5PF.txt")).Data()));
    jetCorr->AddCorrectionFromFile(std::string((gSystem->Getenv("CMSSW_BASE") + TString("/src/MitPhysics/data/Summer13_V1_MC_L3Absolute_AK5PF.txt")).Data()));
  }
  jetCorr->SetInputName(pubJet->GetOutputName());
  jetCorr->SetCorrectedName("CorrectedJets");

  JetIDMod *jetID = new JetIDMod;
  jetID->SetInputName(jetCorr->GetOutputName());
  jetID->SetPtCut(30.0);
  jetID->SetEtaMaxCut(4.7);
  jetID->SetJetEEMFractionMinCut(0.00);
  jetID->SetOutputName("GoodJets");
  jetID->SetApplyBetaCut(kFALSE);
  jetID->SetApplyMVACut(kTRUE);

  JetCleaningMod *jetCleaning = new JetCleaningMod;
  jetCleaning->SetCleanElectronsName(electronCleaning->GetOutputName());
  jetCleaning->SetCleanMuonsName(muonIdIso->GetOutputName()); // clean up isolated muons (instead of the loose ones)
  jetCleaning->SetCleanPhotonsName(photonCleaningMod->GetOutputName());
  jetCleaning->SetApplyPhotonRemoval(kTRUE);
  jetCleaning->SetGoodJetsName(jetID->GetOutputName());
  jetCleaning->SetCleanJetsName("CleanJets");

  MetCorrectionMod *metCorrT0T1Shift = new MetCorrectionMod;
  metCorrT0T1Shift->SetInputName("PFMet");
  metCorrT0T1Shift->SetJetsName(pubJet->GetOutputName());
  metCorrT0T1Shift->SetCorrectedJetsName(jetCorr->GetOutputName());
  metCorrT0T1Shift->SetCorrectedName("PFMetT0T1Shift");
  metCorrT0T1Shift->ApplyType0(kTRUE);
  metCorrT0T1Shift->ApplyType1(kTRUE);
  metCorrT0T1Shift->ApplyShift(kTRUE);
  metCorrT0T1Shift->IsData(isData);
  metCorrT0T1Shift->SetPrint(kFALSE);

  //------------------------------------------------------------------------------------------------
  // select events with jet+MET
  //------------------------------------------------------------------------------------------------

  // VBF
  float minLeadingJetEt = 40;
  float maxJetEta = 4.7;
  float minMet = 110;

  // Monojet
//   float minLeadingJetEt = 40;
//   float maxJetEta = 4.5;
//   float minMet = 200;

  MonoJetAnalysisMod         *jetplusmet = new MonoJetAnalysisMod("MonoJetSelector");
  jetplusmet->SetInputMetName(metCorrT0T1Shift->GetOutputName()); //corrected met
  jetplusmet->SetMetFromBranch(kFALSE);
  jetplusmet->SetJetsName(jetCleaning->GetOutputName()); //identified jets
  jetplusmet->SetJetsFromBranch(kFALSE);
  jetplusmet->SetElectronsName(electronCleaning->GetOutputName());
  jetplusmet->SetElectronsFromBranch(kFALSE);
  jetplusmet->SetMuonsName(muonId->GetOutputName());
  jetplusmet->SetMuonsFromBranch(kFALSE);
  jetplusmet->SetTausName(pftauCleaningMod->GetOutputName());
  jetplusmet->SetTausFromBranch(kFALSE);
  jetplusmet->SetLeptonsName(merger->GetOutputName());
  jetplusmet->SetMinNumJets(1);
  jetplusmet->SetMinNumLeptons(0);
  jetplusmet->SetMinChargedHadronFrac(0.2);
  jetplusmet->SetMaxNeutralHadronFrac(0.7);
  jetplusmet->SetMaxNeutralEmFrac(0.7);
  jetplusmet->SetMinJetEt(minLeadingJetEt); // 40
  jetplusmet->SetMaxJetEta(maxJetEta); //4.7, FIXME: add cut for 2nd jet offline!
  jetplusmet->SetMinMetEt(minMet); // 110, too low?

  MonoJetAnalysisMod         *dilepton = new MonoJetAnalysisMod("MonoJetSelector_dilepton");
  dilepton->SetInputMetName(metCorrT0T1Shift->GetOutputName()); //corrected met
  dilepton->SetMetFromBranch(kFALSE);
  dilepton->SetJetsName(jetCleaning->GetOutputName()); //identified jets
  dilepton->SetJetsFromBranch(kFALSE);
  dilepton->SetElectronsName(electronCleaning->GetOutputName());
  dilepton->SetElectronsFromBranch(kFALSE);
  dilepton->SetMuonsName(muonId->GetOutputName());
  dilepton->SetMuonsFromBranch(kFALSE);
  dilepton->SetTausName(pftauCleaningMod->GetOutputName());
  dilepton->SetTausFromBranch(kFALSE);
  dilepton->SetLeptonsName(merger->GetOutputName());
  dilepton->SetMinNumJets(1);
  dilepton->SetMinNumLeptons(2);
  dilepton->SetMinChargedHadronFrac(0.2);
  dilepton->SetMaxNeutralHadronFrac(0.7);
  dilepton->SetMaxNeutralEmFrac(0.7);
  dilepton->SetMinJetEt(minLeadingJetEt);
  dilepton->SetMaxJetEta(maxJetEta);
  dilepton->SetMinMetEt(0);

  MonoJetAnalysisMod         *wlnu = new MonoJetAnalysisMod("MonoJetSelector_wlnu");
  wlnu->SetInputMetName(metCorrT0T1Shift->GetOutputName()); //corrected met
  wlnu->SetMetFromBranch(kFALSE);
  wlnu->SetJetsName(jetCleaning->GetOutputName()); //identified jets
  wlnu->SetJetsFromBranch(kFALSE);
  wlnu->SetElectronsName(electronCleaning->GetOutputName());
  wlnu->SetElectronsFromBranch(kFALSE);
  wlnu->SetMuonsName(muonId->GetOutputName());
  wlnu->SetMuonsFromBranch(kFALSE);
  wlnu->SetTausName(pftauCleaningMod->GetOutputName());
  wlnu->SetTausFromBranch(kFALSE);
  wlnu->SetLeptonsName(merger->GetOutputName());
  wlnu->SetMinNumJets(1);
  wlnu->SetMinNumLeptons(1);
  wlnu->SetMinChargedHadronFrac(0.2);
  wlnu->SetMaxNeutralHadronFrac(0.7);
  wlnu->SetMaxNeutralEmFrac(0.7);
  wlnu->SetMinJetEt(minLeadingJetEt);
  wlnu->SetMaxJetEta(maxJetEta);
  wlnu->SetMinMetEt(0);

  MonoJetTreeWriter *jetplusmettree = new MonoJetTreeWriter("MonoJetTreeWriter");
  jetplusmettree->SetTriggerObjectsName(hltModP->GetOutputName());
  jetplusmettree->SetMetName(metCorrT0T1Shift->GetOutputName()); //corrected met
  jetplusmettree->SetMetFromBranch(kFALSE);
  jetplusmettree->SetPhotonsFromBranch(kFALSE);
  jetplusmettree->SetPhotonsName(photonCleaningMod->GetOutputName());
  jetplusmettree->SetElectronsFromBranch(kFALSE);
  jetplusmettree->SetElectronsName(electronCleaning->GetOutputName());
  jetplusmettree->SetMuonsFromBranch(kFALSE);
  jetplusmettree->SetMuonsName(muonId->GetOutputName());
  jetplusmettree->SetTausFromBranch(kFALSE);
  jetplusmettree->SetTausName(pftauCleaningMod->GetOutputName());
  jetplusmettree->SetJetsFromBranch(kFALSE);
  jetplusmettree->SetJetsName(jetCleaning->GetOutputName());
  jetplusmettree->SetRawJetsName(pubJet->GetOutputName());
  jetplusmettree->SetPVFromBranch(kFALSE);
  jetplusmettree->SetPVName(goodPVFilterMod->GetOutputName());
  jetplusmettree->SetLeptonsName(merger->GetOutputName());
  jetplusmettree->SetIsData(isData);
  jetplusmettree->SetProcessID(0);
  jetplusmettree->SetFillNtupleType(0);

  MonoJetTreeWriter *dileptontree = new MonoJetTreeWriter("MonoJetTreeWriter_dilepton");
  dileptontree->SetTriggerObjectsName(hltModP->GetOutputName());
  dileptontree->SetMetName(metCorrT0T1Shift->GetOutputName()); //corrected met
  dileptontree->SetMetFromBranch(kFALSE);
  dileptontree->SetPhotonsFromBranch(kFALSE);
  dileptontree->SetPhotonsName(photonCleaningMod->GetOutputName());
  dileptontree->SetElectronsFromBranch(kFALSE);
  dileptontree->SetElectronsName(electronCleaning->GetOutputName());
  dileptontree->SetMuonsFromBranch(kFALSE);
  dileptontree->SetMuonsName(muonId->GetOutputName());
  dileptontree->SetTausFromBranch(kFALSE);
  dileptontree->SetTausName(pftauCleaningMod->GetOutputName());
  dileptontree->SetJetsFromBranch(kFALSE);
  dileptontree->SetJetsName(jetCleaning->GetOutputName());
  dileptontree->SetRawJetsName(pubJet->GetOutputName());
  dileptontree->SetPVFromBranch(kFALSE);
  dileptontree->SetPVName(goodPVFilterMod->GetOutputName());
  dileptontree->SetLeptonsName(merger->GetOutputName());
  dileptontree->SetIsData(isData);
  dileptontree->SetProcessID(0);
  dileptontree->SetFillNtupleType(1);

  MonoJetTreeWriter *wlnutree = new MonoJetTreeWriter("MonoJetTreeWriter_wlnu");
  wlnutree->SetTriggerObjectsName(hltModP->GetOutputName());
  wlnutree->SetMetName(metCorrT0T1Shift->GetOutputName()); //corrected met
  wlnutree->SetMetFromBranch(kFALSE);
  wlnutree->SetPhotonsFromBranch(kFALSE);
  wlnutree->SetPhotonsName(photonCleaningMod->GetOutputName());
  wlnutree->SetElectronsFromBranch(kFALSE);
  wlnutree->SetElectronsName(electronCleaning->GetOutputName());
  wlnutree->SetMuonsFromBranch(kFALSE);
  wlnutree->SetMuonsName(muonId->GetOutputName());
  wlnutree->SetTausFromBranch(kFALSE);
  wlnutree->SetTausName(pftauCleaningMod->GetOutputName());
  wlnutree->SetJetsFromBranch(kFALSE);
  wlnutree->SetJetsName(jetCleaning->GetOutputName());
  wlnutree->SetRawJetsName(pubJet->GetOutputName());
  wlnutree->SetPVFromBranch(kFALSE);
  wlnutree->SetPVName(goodPVFilterMod->GetOutputName());
  wlnutree->SetLeptonsName(merger->GetOutputName());
  wlnutree->SetIsData(isData);
  wlnutree->SetProcessID(0);
  wlnutree->SetFillNtupleType(2);

  //------------------------------------------------------------------------------------------------
  // making analysis chain
  //------------------------------------------------------------------------------------------------
  // this is how it always starts
  runLumiSel->Add(generatorMod);
  generatorMod->Add(goodPVFilterMod);
  goodPVFilterMod->Add(hltModP);
  // photon regression
  hltModP->Add(photreg);
  // simple object id modules
  photreg          ->Add(SepPUMod);
  SepPUMod         ->Add(muonId);
  muonId           ->Add(eleIdMod);
  eleIdMod         ->Add(electronCleaning);
  electronCleaning ->Add(merger);
  merger           ->Add(photonIDMod);
  photonIDMod      ->Add(photonCleaningMod);
  photonCleaningMod->Add(pftauIDMod);
  pftauIDMod       ->Add(pftauCleaningMod);
  pftauCleaningMod ->Add(pubJet);
  pubJet           ->Add(jetCorr);
  jetCorr          ->Add(metCorrT0T1Shift);
  metCorrT0T1Shift ->Add(jetID);
  jetID            ->Add(muonIdIso);

  // Jet+met selection
  muonIdIso        ->Add(jetCleaning);
  jetCleaning      ->Add(jetplusmet);
  jetplusmet       ->Add(jetplusmettree);

  // Dilepton selection
  jetCleaning      ->Add(dilepton);
  dilepton         ->Add(dileptontree);

  // Wlnu selection
  jetCleaning     ->Add(wlnu);
  wlnu            ->Add(wlnutree);

  //------------------------------------------------------------------------------------------------
  // setup analysis
  //------------------------------------------------------------------------------------------------
  Analysis *ana = new Analysis;
  ana->SetUseCacher(cacher);
  ana->SetUseHLT(kTRUE);
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
  Catalog *c = new Catalog(cataDir.Data());
  Dataset *d = NULL;
  if (TString(skim).CompareTo("noskim") == 0)
    d = c->FindDataset(bookstr,dataset,fileset,local);
  else
    d = c->FindDataset(bookstr,skimdataset.Data(),fileset,local);
  ana->AddDataset(d);
  //ana->AddFile("/mnt/hadoop/cms/store/user/paus/filefi/032/r12a-met-j22-v1/C4AC0AB8-BA82-E211-B238-003048678FF4.root");

  //------------------------------------------------------------------------------------------------
  // organize output
  //------------------------------------------------------------------------------------------------
  TString rootFile = TString(outputName);
  rootFile += TString("_") + TString(dataset) + TString("_") + TString(skim);
  if (TString(fileset) != TString(""))
    rootFile += TString("_") + TString(fileset);
  rootFile += TString(".root");
  ana->SetOutputName(rootFile.Data());
  ana->SetCacheSize(0);

  //------------------------------------------------------------------------------------------------
  // Say what we are doing
  //------------------------------------------------------------------------------------------------
  printf("\n==== PARAMETER SUMMARY FOR THIS JOB ====\n");
  printf("\n JSON file: %s\n",jsonFile.Data());
  printf("\n Rely on Catalog: %s\n",cataDir.Data());
  printf("  -> Book: %s  Dataset: %s  Skim: %s  Fileset: %s <-\n",book,dataset,skim,fileset);
  printf("\n Root output: %s\n\n",rootFile.Data());
  printf("\n========================================\n");

  //------------------------------------------------------------------------------------------------
  // run the analysis after successful initialisation
  //------------------------------------------------------------------------------------------------
  ana->Run(kFALSE);

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


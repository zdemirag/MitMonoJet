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
#include "MitAna/DataTree/interface/Names.h"
#include "MitPhysics/Init/interface/ModNames.h"
#include "MitPhysics/Utils/interface/RhoUtilities.h"
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
#include "MitMonoJet/Mods/interface/BoostedVTreeWriter.h"
#endif

//--------------------------------------------------------------------------------------------------
void runBoostedV(const char *fileset    = "0000",
		 const char *skim       = "noskim",
		 const char *dataset    = "s12-ww-v7a", 
		 const char *book       = "t2mit/filefi/032",
		 const char *catalogDir = "/home/cmsprod/catalog",
		 const char *outputName = "boostedv",
		 int         nEvents    = 1000)
{
  //------------------------------------------------------------------------------------------------
  // some parameters get passed through the environment
  //------------------------------------------------------------------------------------------------
  char json[1024], overlap[1024];
  float overlapCut = -1;
  TString mitData;
  TString cataDir("./catalog");
  Long_t *id=0,*size=0,*flags=0,*mt=0;

  printf(" Try local catalog first: %s\n",cataDir.Data());
  if (gSystem->GetPathInfo(cataDir.Data(),id,size,flags,mt) != 0) {
    cataDir = TString(catalogDir);
    if (gSystem->GetPathInfo(cataDir.Data(),id,size,flags,mt) != 0) {
      printf(" Requested local (./catalog) and specified catalog do not exist. EXIT!\n");
      return;
    }
  }
  else {
    printf(" Local catalog exists: %s using this one.\n",cataDir.Data()); 
  }
  if (gSystem->Getenv("MIT_DATA"))
    mitData = gSystem->Getenv("MIT_DATA");
  else {
    printf(" MIT_DATA was not properly defined. EXIT!\n");
    return;
  } 
  if (gSystem->Getenv("MIT_PROD_JSON"))
    sprintf(json,   "%s",gSystem->Getenv("MIT_PROD_JSON"));
  else {
    printf(" JSON file was not properly defined. EXIT!\n");
    return;
  }
  TString jsonFile = TString("/home/cmsprod/cms/json/") + TString(json);
  Bool_t  isData   = ((jsonFile.CompareTo("/home/cmsprod/cms/json/~") != 0));
  
  if (gSystem->Getenv("MIT_PROD_OVERLAP")) {
    sprintf(overlap,"%s",gSystem->Getenv("MIT_PROD_OVERLAP"));
    if (EOF == sscanf(overlap,"%f",&overlapCut)) {
      printf(" Overlap was not properly defined. EXIT!\n");
      return;
    }
  }
  else {
    printf(" OVERLAP file was not properly defined. EXIT!\n");
    return;
  } 

  printf("\n Initialization worked. Data?: %d\n\n",isData);

  //------------------------------------------------------------------------------------------------
  // some global setups
  //------------------------------------------------------------------------------------------------
  using namespace mithep;
  gDebugMask  = Debug::kGeneral;
  gDebugLevel = 3;

  //------------------------------------------------------------------------------------------------
  // set up information for master module
  //------------------------------------------------------------------------------------------------
  RunLumiSelectionMod *runLumiSel = new RunLumiSelectionMod;
  runLumiSel->SetAcceptMC(kTRUE);                          // Monte Carlo events are always accepted
  
  // only select on run- and lumisection numbers when valid json file present
  if ((jsonFile.CompareTo("/home/cmsprod/cms/json/~") != 0) &&
      (jsonFile.CompareTo("/home/cmsprod/cms/json/-") != 0)   ) {
    runLumiSel->AddJSONFile(jsonFile.Data());
  }
  if ((jsonFile.CompareTo("/home/cmsprod/cms/json/-") == 0)   ) {
    printf("\n WARNING -- Looking at data without JSON file: always accept.\n\n");
    runLumiSel->SetAbortIfNotAccepted(kFALSE);   // accept all events if there is no valid JSON file
  }
  printf("\n Run lumi worked. \n\n");

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
  Catalog *c = new Catalog(cataDir.Data());
  TString skimdataset = TString(dataset)+TString("/") +TString(skim);
  Dataset *d = NULL;
  TString bookstr = book;
  if (TString(skim).CompareTo("noskim") == 0)
    d = c->FindDataset(bookstr,dataset,fileset,true);
  else 
    d = c->FindDataset(bookstr,skimdataset.Data(),fileset,true);
  ana->AddDataset(d);
  //ana->AddFile("/mnt/hadoop/cms/store/user/paus/filefi/032/s12-wjets-ptw100-v7a/FCEC7DAE-A367-E211-A486-0024E87687BE.root");

  //------------------------------------------------------------------------------------------------
  // organize output
  //------------------------------------------------------------------------------------------------
  TString rootFile = TString(outputName);
  rootFile += TString("_") + TString(dataset) + TString("_") + TString(skim);
  if (TString(fileset) != TString(""))
    rootFile += TString("_") + TString(fileset);
  //TString ntupleFile = rootFile + TString("_ntuple.root");
  rootFile += TString(".root");
  ana->SetOutputName(rootFile.Data());
  ana->SetCacheSize(0);

  //------------------------------------------------------------------------------------------------
  // HLT information
  //------------------------------------------------------------------------------------------------
  HLTMod *hltModP = new HLTMod("HLTModP");
  hltModP->SetBitsName("HLTBits");
  hltModP->SetTrigObjsName("HltObjsMonoJet");
  hltModP->SetAbortIfNotAccepted(isData);
  hltModP->SetPrintTable(kFALSE);

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

  //------------------------------------------------------------------------------------------------
  // split pfcandidates to PFPU and PFnoPU
  //------------------------------------------------------------------------------------------------
  SeparatePileUpMod* sepPuMod = new SeparatePileUpMod;
  //sepPuMod->SetUseAllVerteces(kFALSE);
  //sepPuMod->SetVertexName("OutVtxCiC");
  sepPuMod->SetPFNoPileUpName("pfnopileupcands");
  sepPuMod->SetPFPileUpName("pfpileupcands");
  sepPuMod->SetCheckClosestZVertex(kFALSE);

  //------------------------------------------------------------------------------------------------
  // select events with a good primary vertex
  //------------------------------------------------------------------------------------------------
  GoodPVFilterMod *goodPVFilterMod = new GoodPVFilterMod;
  goodPVFilterMod->SetMinVertexNTracks(0);
  goodPVFilterMod->SetMinNDof(4.0);
  goodPVFilterMod->SetMaxAbsZ(24.0);
  goodPVFilterMod->SetMaxRho(2.0);
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
  eleIdMod->SetIDType("VBTFWorkingPoint95Id");
  eleIdMod->SetIsoType("PFIso");
  eleIdMod->SetApplyConversionFilterType1(kTRUE);
  eleIdMod->SetApplyConversionFilterType2(kFALSE);
  eleIdMod->SetChargeFilter(kFALSE);
  eleIdMod->SetApplyD0Cut(kTRUE);
  eleIdMod->SetApplyDZCut(kTRUE);
  eleIdMod->SetWhichVertex(0);
  eleIdMod->SetNExpectedHitsInnerCut(0);
  eleIdMod->SetGoodElectronsName("GoodElectronsBS");
  eleIdMod->SetRhoType(RhoUtilities::CMS_RHO_RHOKT6PFJETS);

  MuonIDMod* muonIdGammaGamma = new MuonIDMod;
  // base kinematics
  muonIdGammaGamma->SetPtMin(10.);
  muonIdGammaGamma->SetEtaCut(2.4);
  // base ID
  muonIdGammaGamma->SetIDType("NoId");
  muonIdGammaGamma->SetClassType("GlobalorTracker");
  muonIdGammaGamma->SetWhichVertex(0); // this is a 'hack'.. but hopefully good enough...
  muonIdGammaGamma->SetD0Cut(0.02);
  muonIdGammaGamma->SetDZCut(0.5);
  muonIdGammaGamma->SetIsoType("PFIsoBetaPUCorrected"); //h
  muonIdGammaGamma->SetPFIsoCut(0.2); //h
  muonIdGammaGamma->SetOutputName("HggLeptonTagMuons");
  muonIdGammaGamma->SetPFNoPileUpName("pfnopileupcands");
  muonIdGammaGamma->SetPFPileUpName("pfpileupcands");
  muonIdGammaGamma->SetPVName(Names::gkPVBeamSpotBrn);

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

  MuonIDMod *muonId = muonIdWW;
  //MuonIDMod *muonId = muonIdGammaGamma;

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
  PhotonMvaMod *photonReg = new PhotonMvaMod;
  photonReg->SetRegressionVersion(3);
  photonReg->SetRegressionWeights((mitData+TString("/gbrv3ph_52x.root")).Data());
  photonReg->SetOutputName("GoodPhotonsRegr");
  photonReg->SetApplyShowerRescaling(kTRUE);
  photonReg->SetMinNumPhotons(0);
  photonReg->SetIsData(isData);

  PhotonIDMod *photonIdMod = new PhotonIDMod;
  photonIdMod->SetPtMin(0.0);
  photonIdMod->SetOutputName("GoodPhotons");
  photonIdMod->SetIDType("BaseLineCiCPFNoPresel");
  photonIdMod->SetIsoType("NoIso");
  photonIdMod->SetApplyElectronVeto(kTRUE);
  photonIdMod->SetApplyPixelSeed(kTRUE);
  photonIdMod->SetApplyConversionId(kTRUE);
  photonIdMod->SetApplyFiduciality(kTRUE);
  photonIdMod->SetIsData(isData);
  photonIdMod->SetPhotonsFromBranch(kFALSE);
  photonIdMod->SetInputName(photonReg->GetOutputName());
  //get the photon with regression energy
  photonIdMod->DoMCSmear(kTRUE);
  photonIdMod->DoDataEneCorr(kTRUE);
  //------------------------------------------Energy smear and scale--------------------------------------------------------------
  photonIdMod->SetMCSmearFactors2012HCP(0.0111,0.0111,0.0107,0.0107,0.0155,0.0194,0.0295,0.0276,0.037,0.0371);
  photonIdMod->AddEnCorrPerRun2012HCP(190645,190781,0.9964,0.9964,1.0020,1.0020,0.9893,1.0028,0.9871,0.9937,0.9839,0.9958);
  photonIdMod->AddEnCorrPerRun2012HCP(190782,191042,1.0024,1.0024,1.0079,1.0079,0.9923,1.0058,0.9911,0.9977,0.9886,1.0005);
  photonIdMod->AddEnCorrPerRun2012HCP(191043,193555,0.9935,0.9935,0.9991,0.9991,0.9861,0.9997,0.9894,0.9960,0.9864,0.9982);
  photonIdMod->AddEnCorrPerRun2012HCP(193556,194150,0.9920,0.9920,0.9976,0.9976,0.9814,0.9951,0.9896,0.9962,0.9872,0.9990);
  photonIdMod->AddEnCorrPerRun2012HCP(194151,194532,0.9925,0.9925,0.9981,0.9981,0.9826,0.9963,0.9914,0.9980,0.9874,0.9993);
  photonIdMod->AddEnCorrPerRun2012HCP(194533,195113,0.9927,0.9927,0.9983,0.9983,0.9844,0.9981,0.9934,0.9999,0.9878,0.9996);
  photonIdMod->AddEnCorrPerRun2012HCP(195114,195915,0.9929,0.9929,0.9984,0.9984,0.9838,0.9974,0.9942,1.0007,0.9878,0.9997);
  photonIdMod->AddEnCorrPerRun2012HCP(195916,198115,0.9919,0.9919,0.9975,0.9975,0.9827,0.9964,0.9952,1.0017,0.9869,0.9987);
  photonIdMod->AddEnCorrPerRun2012HCP(198116,199803,0.9955,0.9955,1.0011,1.0011,0.9859,0.9995,0.9893,0.9959,0.9923,1.0041);
  photonIdMod->AddEnCorrPerRun2012HCP(199804,200048,0.9967,0.9967,1.0023,1.0023,0.9870,1.0006,0.9893,0.9959,0.9937,1.0055);
  photonIdMod->AddEnCorrPerRun2012HCP(200049,200151,0.9980,0.9980,1.0036,1.0036,0.9877,1.0012,0.9910,0.9976,0.9980,1.0097);
  photonIdMod->AddEnCorrPerRun2012HCP(200152,200490,0.9958,0.9958,1.0013,1.0013,0.9868,1.0004,0.9922,0.9988,0.9948,1.0065);
  photonIdMod->AddEnCorrPerRun2012HCP(200491,200531,0.9979,0.9979,1.0035,1.0035,0.9876,1.0012,0.9915,0.9981,0.9979,1.0096);
  photonIdMod->AddEnCorrPerRun2012HCP(200532,201656,0.9961,0.9961,1.0017,1.0017,0.9860,0.9996,0.9904,0.9970,0.9945,1.0063);
  photonIdMod->AddEnCorrPerRun2012HCP(201657,202305,0.9969,0.9969,1.0025,1.0025,0.9866,1.0002,0.9914,0.9980,0.9999,1.0116);
  photonIdMod->AddEnCorrPerRun2012HCP(202305,203002,0.9982,0.9982,1.0038,1.0038,0.9872,1.0008,0.9934,1.0000,1.0018,1.0135);
  photonIdMod->AddEnCorrPerRun2012HCP(203003,203984,1.0006,1.0006,1.0061,1.0061,0.9880,1.0017,0.9919,0.9988,0.9992,1.0104);
  photonIdMod->AddEnCorrPerRun2012HCP(203985,205085,0.9993,0.9993,1.0048,1.0048,0.9903,1.0040,0.9928,0.9997,0.9987,1.0099);
  photonIdMod->AddEnCorrPerRun2012HCP(205086,205310,1.0004,1.0004,1.0059,1.0059,0.9901,1.0037,0.9987,1.0055,1.0091,1.0202);
  photonIdMod->AddEnCorrPerRun2012HCP(205311,206207,1.0000,1.0000,1.0055,1.0055,0.9891,1.0028,0.9948,1.0017,1.0032,1.0144);
  photonIdMod->AddEnCorrPerRun2012HCP(206208,206483,1.0003,1.0003,1.0058,1.0058,0.9895,1.0032,0.9921,0.9989,1.0056,1.0167);
  photonIdMod->AddEnCorrPerRun2012HCP(206484,206597,1.0005,1.0005,1.0060,1.0060,0.9895,1.0032,0.9968,1.0036,1.0046,1.0158);
  photonIdMod->AddEnCorrPerRun2012HCP(206598,206896,1.0006,1.0006,1.0061,1.0061,0.9881,1.0017,0.9913,0.9982,1.0050,1.0162);
  photonIdMod->AddEnCorrPerRun2012HCP(206897,207220,1.0006,1.0006,1.0061,1.0061,0.9884,1.0021,0.9909,0.9978,1.0053,1.0165);
  photonIdMod->AddEnCorrPerRun2012HCP(207221,208686,1.0006,1.0006,1.0061,1.0061,0.9894,1.0030,0.9951,1.0020,1.0060,1.0172);
  //---------------------------------shower shape scale--------------------------------------------------------------------------------
  photonIdMod->SetDoShowerShapeScaling(kTRUE);
  photonIdMod->SetShowerShapeType("2012ShowerShape");
  photonIdMod->Set2012HCP(kTRUE);

  PFTauIDMod *pftauIdMod = new PFTauIDMod;
  pftauIdMod->SetPFTausName("HPSTaus");
  pftauIdMod->SetIsLooseId(kFALSE);

  PhotonCleaningMod *photonCleaningMod = new PhotonCleaningMod;
  photonCleaningMod->SetCleanElectronsName(electronCleaning->GetOutputName());
  photonCleaningMod->SetGoodPhotonsName(photonIdMod->GetOutputName());
  photonCleaningMod->SetCleanPhotonsName("CleanPhotons");

  PFTauCleaningMod *pftauCleaningMod = new PFTauCleaningMod;
  pftauCleaningMod->SetGoodPFTausName(pftauIdMod->GetGoodPFTausName());
  pftauCleaningMod->SetCleanMuonsName(muonId->GetOutputName());

  PublisherMod<PFJet,Jet> *pubJet = new PublisherMod<PFJet,Jet>("JetPub");
  pubJet->SetInputName("AKt5PFJets");
  pubJet->SetOutputName("PubAKt5PFJets");

  JetCorrectionMod *jetCorr = new JetCorrectionMod;
  if (isData){ 
    jetCorr->AddCorrectionFromFile((mitData+TString("/Summer13_V1_DATA_L1FastJet_AK5PF.txt")).Data()); 
    jetCorr->AddCorrectionFromFile((mitData+TString("/Summer13_V1_DATA_L2Relative_AK5PF.txt")).Data()); 
    jetCorr->AddCorrectionFromFile((mitData+TString("/Summer13_V1_DATA_L3Absolute_AK5PF.txt")).Data()); 
    jetCorr->AddCorrectionFromFile((mitData+TString("/Summer13_V1_DATA_L2L3Residual_AK5PF.txt")).Data());
  }                                                                                      
  else {                                                                                 
    jetCorr->AddCorrectionFromFile((mitData+TString("/Summer13_V1_MC_L1FastJet_AK5PF.txt")).Data()); 
    jetCorr->AddCorrectionFromFile((mitData+TString("/Summer13_V1_MC_L2Relative_AK5PF.txt")).Data()); 
    jetCorr->AddCorrectionFromFile((mitData+TString("/Summer13_V1_MC_L3Absolute_AK5PF.txt")).Data()); 
  }
  jetCorr->SetInputName(pubJet->GetOutputName());
  jetCorr->SetCorrectedName("CorrectedJets");    
        
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

  JetIDMod *jetId = new JetIDMod;
  jetId->SetInputName(jetCorr->GetOutputName());
  jetId->SetPtCut(30.0);
  jetId->SetEtaMaxCut(4.7);
  jetId->SetJetEEMFractionMinCut(0.00);
  jetId->SetOutputName("GoodJets");
  jetId->SetApplyBetaCut(kFALSE);
  jetId->SetApplyMVACut(kTRUE);

  JetCleaningMod *jetCleaning = new JetCleaningMod;
  jetCleaning->SetCleanElectronsName(electronCleaning->GetOutputName());
  jetCleaning->SetCleanMuonsName(muonId->GetOutputName());
  jetCleaning->SetCleanPhotonsName(photonCleaningMod->GetOutputName());
  jetCleaning->SetApplyPhotonRemoval(kTRUE);
  jetCleaning->SetGoodJetsName(jetId->GetOutputName());
  jetCleaning->SetCleanJetsName("CleanJets");
 
  //------------------------------------------------------------------------------------------------
  // select events with a given jet substructure
  //------------------------------------------------------------------------------------------------
  BoostedVTreeWriter *boostedVMod = new BoostedVTreeWriter;
  boostedVMod->SetIsData(isData);
  boostedVMod->SetTriggerObjsName(hltModP->GetOutputName());
  boostedVMod->SetJetsName(jetCleaning->GetOutputName());
  boostedVMod->SetJetsFromBranch(kFALSE);
  boostedVMod->SetPhotonsName(photonCleaningMod->GetOutputName());
  boostedVMod->SetPhotonsFromBranch(kFALSE);
  boostedVMod->SetPFTausName(pftauCleaningMod->GetOutputName());
  boostedVMod->SetPFTausFromBranch(kFALSE);
  boostedVMod->SetLeptonsName(merger->GetOutputName());
  boostedVMod->SetPruning(1);
  //boostedVMod->SetOutputName(ntupleFile.Data());

  //------------------------------------------------------------------------------------------------
  // making analysis chain
  //------------------------------------------------------------------------------------------------
  runLumiSel       ->Add(goodPVFilterMod);
  goodPVFilterMod  ->Add(hltModP);
  hltModP          ->Add(photonReg);
  photonReg        ->Add(sepPuMod);
  sepPuMod         ->Add(muonId);
  muonId           ->Add(eleIdMod);
  eleIdMod         ->Add(electronCleaning);
  electronCleaning ->Add(merger);
  merger           ->Add(photonIdMod);
  photonIdMod      ->Add(photonCleaningMod);
  photonCleaningMod->Add(pftauIdMod);
  pftauIdMod       ->Add(pftauCleaningMod);
  pftauCleaningMod ->Add(pubJet);
  pubJet           ->Add(jetCorr);
  jetCorr          ->Add(metCorrT0T1Shift);
  metCorrT0T1Shift ->Add(jetId);
  jetId            ->Add(jetCleaning);
  jetCleaning      ->Add(boostedVMod);
  
  //------------------------------------------------------------------------------------------------
  // Say what we are doing
  //------------------------------------------------------------------------------------------------
  printf("\n==== PARAMETER SUMMARY FOR THIS JOB ====\n");
  printf("\n JSON file: %s\n  and overlap cut: %f (%s)\n",jsonFile.Data(),overlapCut,overlap);
  printf("\n Rely on Catalog: %s\n",cataDir.Data());
  printf("  -> Book: %s  Dataset: %s  Skim: %s  Fileset: %s <-\n",book,dataset,skim,fileset);
  printf("\n Root output:   %s\n",rootFile.Data());  
  //printf("\n Ntuple output: %s\n\n",ntupleFile.Data());  
  printf("\n========================================\n");

  //------------------------------------------------------------------------------------------------
  // run the analysis after successful initialisation
  //------------------------------------------------------------------------------------------------
  ana->Run(!gROOT->IsBatch());

  return;
}

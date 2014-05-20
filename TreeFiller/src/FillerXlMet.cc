#include <TSystem.h> 

#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"

#include "MitAna/DataTree/interface/MuonFwd.h"
#include "MitAna/DataTree/interface/ElectronFwd.h"
#include "MitAna/DataTree/interface/VertexCol.h"
#include "MitAna/DataTree/interface/PFCandidateCol.h"
#include "MitAna/DataTree/interface/PFMet.h"
#include "MitAna/DataTree/interface/PFJetCol.h"
#include "MitAna/DataTree/interface/PFTauCol.h"
#include "MitPhysics/Init/interface/ModNames.h"
#include "MitPhysics/Utils/interface/MuonTools.h"
#include "MitMonoJet/TreeFiller/interface/FillerXlMet.h"
#include "MitMonoJet/DataTree/interface/XlMet.h"

using namespace mithep;

ClassImp(mithep::FillerXlMet)

//--------------------------------------------------------------------------------------------------
FillerXlMet::FillerXlMet(const char *name, const char *title) :
  BaseMod (name,title),
  fIsData (kTRUE),
  fJetsName (Names::gkPFJetBrn),
  fMuonsName (Names::gkMuonBrn),
  fElectronsName (Names::gkElectronBrn),
  fTausName (Names::gkPFTauBrn),
  fPFCandidatesName (Names::gkPFCandidatesBrn),
  fPVName (Names::gkPVBeamSpotBrn),
  fPileUpDenName (Names::gkPileupEnergyDensityBrn),
  fRawMetName ("PFMet"),
  fXlMetName ("XlMetMVA"),
  fJetsFromBranch (kTRUE),
  fMuonsFromBranch (kTRUE),
  fElectronsFromBranch (kTRUE),
  fTausFromBranch (kTRUE),
  fPFCandidatesFromBranch (kTRUE),
  fPVFromBranch (kTRUE),
  fPublishOutput (kTRUE),
  fJetCorrector (0)
{
  // Constructor.
}

FillerXlMet::~FillerXlMet()
{
  // Destructor
  if (fXlMet)
    delete fXlMet;
}

//--------------------------------------------------------------------------------------------------
void FillerXlMet::Process()
{
  // make sure the out collections are empty before starting
  fXlMet->Delete();  
  
  // Load the branches we want to work with
  LoadEventObject(fJetsName, fJets, fJetsFromBranch);
  LoadEventObject(fMuonsName, fMuons, fMuonsFromBranch);
  LoadEventObject(fElectronsName, fElectrons, fElectronsFromBranch);
  LoadEventObject(fTausName, fPFTaus, fTausFromBranch);
  LoadEventObject(fPFCandidatesName,  fPFCandidates,  fPFCandidatesFromBranch);
  LoadEventObject(fPVName, fPV, fPVFromBranch);
  LoadEventObject(fRawMetName, fRawMet, true);
  LoadEventObject(fPileUpDenName, fPileUpDen, true);

  // Convert the input collection into PFJets
  PFJetOArr *fPFJets = new PFJetOArr;
  for(UInt_t i=0; i<fJets->GetEntries(); i++) {
    const PFJet *pfJet = dynamic_cast<const PFJet*>(fJets->At(i));
    fPFJets->Add(pfJet);
  }

  // Compute the MVA met for this event
  Met mvaMet = fMVAMet->GetMet(fMuons,fElectrons,fPFTaus,fPFCandidates,
                               fPFJets,0,fPV,fRawMet,fJetCorrector,fPileUpDen);
  TMatrixD* MVACov = fMVAMet->GetMetCovariance();

  // Prepare and store in an array a new XlMet 
  XlMet *extMet = fXlMet->Allocate();
  new (extMet) XlMet(mvaMet.Mex(),mvaMet.Mey());

  // Store the covariance matrix in the new object
  extMet->SetCovMatrix(*MVACov);
  
  return;
}

//--------------------------------------------------------------------------------------------------
void FillerXlMet::SlaveBegin()
{
  // Run startup code on the computer (slave) doing the actual analysis.
  ReqEventObject(fJetsName, fJets, fJetsFromBranch);
  ReqEventObject(fMuonsName, fMuons, fMuonsFromBranch);
  ReqEventObject(fElectronsName, fElectrons, fElectronsFromBranch);
  ReqEventObject(fTausName, fPFTaus, fTausFromBranch);
  ReqEventObject(fPFCandidatesName,  fPFCandidates,  fPFCandidatesFromBranch);
  ReqEventObject(fPVName, fPV, fPVFromBranch);
  ReqEventObject(fRawMetName, fRawMet, true);
  ReqEventObject(fPileUpDenName, fPileUpDen, true);

  // Setup everything needed for JetCorrector
  if (fIsData) {
    fCorrectionFiles.push_back(std::string((gSystem->Getenv("CMSSW_BASE") + TString("/src/MitPhysics/data/Summer13_V1_DATA_L1FastJet_AK5PF.txt")).Data()));
    fCorrectionFiles.push_back(std::string((gSystem->Getenv("CMSSW_BASE") + TString("/src/MitPhysics/data/Summer13_V1_DATA_L2Relative_AK5PF.txt")).Data()));
    fCorrectionFiles.push_back(std::string((gSystem->Getenv("CMSSW_BASE") + TString("/src/MitPhysics/data/Summer13_V1_DATA_L3Absolute_AK5PF.txt")).Data()));
    fCorrectionFiles.push_back(std::string((gSystem->Getenv("CMSSW_BASE") + TString("/src/MitPhysics/data/Summer13_V1_DATA_L2L3Residual_AK5PF.txt")).Data()));
  }
  else {
    fCorrectionFiles.push_back(std::string((gSystem->Getenv("CMSSW_BASE") + TString("/src/MitPhysics/data/Summer13_V1_MC_L1FastJet_AK5PF.txt")).Data()));
    fCorrectionFiles.push_back(std::string((gSystem->Getenv("CMSSW_BASE") + TString("/src/MitPhysics/data/Summer13_V1_MC_L2Relative_AK5PF.txt")).Data()));
    fCorrectionFiles.push_back(std::string((gSystem->Getenv("CMSSW_BASE") + TString("/src/MitPhysics/data/Summer13_V1_MC_L3Absolute_AK5PF.txt")).Data()));
  }

  // Initialize JetCorrectorParameters from files
  std::vector<JetCorrectorParameters> correctionParameters;
  for (std::vector<std::string>::const_iterator it = fCorrectionFiles.begin(); it!=fCorrectionFiles.end(); ++it)
    correctionParameters.push_back(JetCorrectorParameters(*it));

  // Initialize jet corrector class
  fJetCorrector = new FactorizedJetCorrector(correctionParameters);

  // Create a new MVA MET object
  fMVAMet = new MVAMet();  
  fMVAMet->Initialize(
  TString(getenv("CMSSW_BASE")+string("/src/MitPhysics/data/TMVAClassificationCategory_JetID_MET_53X_Dec2012.weights.xml")),
  TString(getenv("CMSSW_BASE")+string("/src/MitPhysics/data/TMVAClassificationCategory_JetID_MET_53X_Dec2012.weights.xml")),
  TString(getenv("CMSSW_BASE")+string("/src/MitPhysics/Utils/python/JetIdParams_cfi.py")),
  TString(getenv("CMSSW_BASE")+string("/src/MitPhysics/data/gbrmet_53_June2013_type1.root")),
  TString(getenv("CMSSW_BASE")+string("/src/MitPhysics/data/gbrmetphi_53_June2013_type1.root")),
  TString(getenv("CMSSW_BASE")+string("/src/MitPhysics/data/gbru1cov_53_Dec2012.root")),
  TString(getenv("CMSSW_BASE")+string("/src/MitPhysics/data/gbru2cov_53_Dec2012.root")),JetIDMVA::k53MET,MVAMet::kUseType1Rho
  );

  // Create the new output collection
  fXlMet = new XlMetArr(16,fXlMetName);
  // Publish collection for further usage in the analysis
  if (fPublishOutput)
    PublishObj(fXlMet);
  
  return;
}

//--------------------------------------------------------------------------------------------------
void FillerXlMet::SlaveTerminate()
{
  delete fMVAMet;
}

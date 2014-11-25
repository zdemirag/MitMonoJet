#include "MitPhysics/Init/interface/ModNames.h"
#include "MitMonoJet/TreeFiller/interface/FillerXsIsoParticles.h"

#include "MitMonoJet/DataTree/interface/XsIsoParticle.h"
#include "MitCommon/DataFormats/interface/Vect4M.h"
#include "MitCommon/DataFormats/interface/Vect3.h"
#include "MitCommon/DataFormats/interface/Types.h"
#include "MitCommon/MathTools/interface/MathUtils.h"

using namespace mithep;

ClassImp(mithep::FillerXsIsoParticles)

//--------------------------------------------------------------------------------------------------
FillerXsIsoParticles::FillerXsIsoParticles(const char *name, const char *title) :
  BaseMod (name,title),
  fIsData (kTRUE),
  fFillXsMuons (kTRUE),
  fFillXsElectrons (kTRUE),
  fFillXsTaus (kTRUE),
  fFillXsPhotons (kTRUE),
  fPublishOutput (kTRUE),
  fMuonsName (Names::gkMuonBrn),
  fMuonsFromBranch (kTRUE),
  fMuons (0),
  fGoodMuonsName (Names::gkMuonBrn),
  fGoodMuonsFromBranch (kFALSE),
  fGoodMuons (0),
  fElectronsName (Names::gkElectronBrn),
  fElectronsFromBranch (kTRUE),
  fElectrons (0),
  fTausName (Names::gkPFTauBrn),
  fTausFromBranch (kTRUE),
  fTaus (0),
  fPhotonsName (Names::gkPhotonBrn),
  fPhotonsFromBranch (kTRUE),
  fPhotons (0),
  fXsMuonsName ("XsMuons"),
  fXsElectronsName ("XsElectrons"),
  fXsTausName ("XsTaus"),
  fXsPhotonsName ("XsPhotons")
{
  // Constructor.
}

FillerXsIsoParticles::~FillerXsIsoParticles()
{
  // Destructor
  if (fXsMuons)
    delete fXsMuons;
  if (fXsElectrons)
    delete fXsElectrons;
  if (fXsTaus)
    delete fXsTaus;
  if (fXsPhotons)
    delete fXsPhotons;
}

//--------------------------------------------------------------------------------------------------
void FillerXsIsoParticles::Process()
{
  // make sure the out collections are empty before starting
  fXsMuons->Delete();  
  fXsElectrons->Delete();  
  fXsTaus->Delete();  
  fXsPhotons->Delete();  
  
  // Load the branches we want to work with
  if (fFillXsMuons) {
    LoadEventObject(fMuonsName,fMuons,fMuonsFromBranch);
    LoadEventObject(fGoodMuonsName,fGoodMuons,fGoodMuonsFromBranch);
  }
  if (fFillXsElectrons) 
    LoadEventObject(fElectronsName,fElectrons,fElectronsFromBranch);
  if (fFillXsTaus) 
    LoadEventObject(fTausName,fTaus,fTausFromBranch);
  if (fFillXsPhotons) 
    LoadEventObject(fPhotonsName,fPhotons,fPhotonsFromBranch);

  // Muon block
  if (fFillXsMuons)
    for (UInt_t i=0; i<fMuons->GetEntries(); ++i) {
      const Particle *muon = dynamic_cast<const Particle*>(fMuons->At(i));
      // Fill the XsMuons collection with the reduced muon object
      FillXsIsoParticle(fXsMuons,muon);            
    }

  // Trim the output collections
  if (fFillXsMuons)
    fXsMuons->Trim();  
  if (fFillXsElectrons) 
    fXsElectrons->Trim();  
  if (fFillXsTaus) 
    fXsTaus->Trim();  
  if (fFillXsPhotons) 
    fXsPhotons->Trim();  
    
  return;
}

//--------------------------------------------------------------------------------------------------
void FillerXsIsoParticles::SlaveBegin()
{
  // Run startup code on the computer (slave) doing the actual analysis. 
  if (fFillXsMuons) {
    ReqEventObject(fMuonsName,fMuons,fMuonsFromBranch);
    ReqEventObject(fGoodMuonsName,fGoodMuons,fGoodMuonsFromBranch);
  }
  if (fFillXsElectrons) 
    ReqEventObject(fElectronsName,fElectrons,fElectronsFromBranch);
  if (fFillXsTaus) 
    ReqEventObject(fTausName,fTaus,fTausFromBranch);
  if (fFillXsPhotons) 
    ReqEventObject(fPhotonsName,fPhotons,fPhotonsFromBranch);

  // Create the new output collection
  fXsMuons = new XsIsoParticleArr(16,fXsMuonsName);
  fXsElectrons = new XsIsoParticleArr(16,fXsElectronsName);
  fXsTaus = new XsIsoParticleArr(16,fXsTausName);
  fXsPhotons = new XsIsoParticleArr(16,fXsPhotonsName);
  // Publish collection for further usage in the analysis
  if (fPublishOutput) {
    PublishObj(fXsMuons);
    PublishObj(fXsElectrons);
    PublishObj(fXsTaus);
    PublishObj(fXsPhotons);
  }
  
  return;
}

//--------------------------------------------------------------------------------------------------
void FillerXsIsoParticles::SlaveTerminate()
{
}

//--------------------------------------------------------------------------------------------------
void FillerXsIsoParticles::FillXsIsoParticle(XsIsoParticleArr *pXsArr, const Particle *pParticle)
{
  // Prepare and store in an array a new XsIsoParticle 
  XsIsoParticle *thisXsIsoParticle = pXsArr->Allocate();
  new (thisXsIsoParticle) XsIsoParticle(*pParticle);
  
  return;
}

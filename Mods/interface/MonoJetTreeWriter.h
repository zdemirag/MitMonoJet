//--------------------------------------------------------------------------------------------------
// $Id: MonoJetTreeWriter.h,v 1.9 2013/09/30 23:08:36 mzanetti Exp $
//
// MonoJetTreeWriter
//
// Authors: L. Di Matteo
//--------------------------------------------------------------------------------------------------

#ifndef MITMONOJET_MODS_MonoJetTreeWriter_H
#define MITMONOJET_MODS_MonoJetTreeWriter_H

#include "MitAna/TreeMod/interface/BaseMod.h" 
#include "MitAna/DataTree/interface/PhotonFwd.h"
#include "MitAna/DataTree/interface/TrackCol.h"
#include "MitAna/DataTree/interface/VertexCol.h"
#include "MitAna/DataTree/interface/BeamSpotCol.h"
#include "MitAna/DataTree/interface/PFCandidateCol.h"
#include "MitAna/DataTree/interface/PileupEnergyDensityCol.h"
#include "MitAna/DataTree/interface/DecayParticleCol.h"
#include "MitAna/DataTree/interface/ElectronCol.h"
#include "MitAna/DataTree/interface/MuonCol.h"
#include "MitAna/DataTree/interface/PFTauCol.h"
#include "MitAna/DataTree/interface/DecayParticleCol.h"
#include "MitAna/DataTree/interface/PileupInfoCol.h"
#include "MitAna/DataTree/interface/MCParticleCol.h"
#include "MitAna/DataTree/interface/MCEventInfo.h"
#include "MitAna/DataTree/interface/SuperClusterCol.h"
#include "MitAna/DataTree/interface/MetCol.h"
#include "MitAna/DataTree/interface/PFMetCol.h"
#include "MitAna/DataTree/interface/JetCol.h"
#include "MitAna/DataTree/interface/PFJetCol.h"
#include "MitAna/DataTree/interface/GenJetCol.h"
#include "MitAna/DataTree/interface/MCEventInfo.h"
#include "MitAna/DataTree/interface/MCParticleFwd.h"
#include "MitPhysics/Utils/interface/PhotonFix.h"
#include "MitPhysics/Utils/interface/PhotonTools.h"
#include "MitPhysics/Utils/interface/MVAMet.h"
#include "MitPhysics/Utils/interface/MVAVBF.h"
#include "MitPhysics/Utils/interface/QGTagger.h"
#include "MitAna/DataTree/interface/EvtSelData.h"

#include "MitPhysics/Utils/interface/VertexTools.h"
#include "MitPhysics/Utils/interface/ElectronIDMVA.h"

#include "MitMonoJet/Core/MitGPTree.h"

class TNtuple;
class TRandom3;

namespace mithep 
{

  class MonoJetTreeWriter : public BaseMod
  {
  public:
    MonoJetTreeWriter(const char *name ="MonoJetTreeWriter", 
		         const char *title="Selecting PhotonPairs");
    
    ~MonoJetTreeWriter();

    // setting all the input Names
    void                SetMetName(const char *n)         { fMetName= n;                 }
    void                SetMetFromBranch(bool b)          { fMetFromBranch = b;          }
    void                SetPhotonsName(const char *n)     { fPhotonsName= n;             }
    void                SetPhotonsFromBranch(bool b)      { fPhotonsFromBranch = b;      }
    void                SetElectronsName(const char *n)   { fElectronsName = n;          }
    void                SetElectronsFromBranch(bool b)    { fElectronsFromBranch = b;    }
    void                SetMuonsName(const char *n)       { fMuonsName = n;              }
    void                SetMuonsFromBranch(bool b)        { fMuonsFromBranch = b;        }
    void                SetTausName(const char *n)        { fTausName = n;               }
    void                SetTausFromBranch(bool b)         { fTausFromBranch = b;         }

    void                SetJetsName(const char *n)        { fJetsName = n;               }
    void                SetJetsFromBranch(bool b)         { fJetsFromBranch = b;         }
    void                SetQGTaggerCHS(bool b)            { fQGTaggerCHS = b;            }
    void                SetLeptonsName(const char *n)     { fLeptonsName = n;            }

    void                SetSuperClustersName(const char *n){ fSuperClustersName = n;     }
    void                SetTracksName(const char *n)      { fTracksName = n;             }
    void                SetPVName(const char *n)          { fPVName = n;                 }
    void                SetPVFromBranch(bool b)           { fPVFromBranch = b;           }
    void                SetPUInfoName(const char *n)      { fPileUpName = n;             }
    void                SetBeamspotName(const char *n)    { fBeamspotName = n;           }

    void                SetIsData (Bool_t b)              { fIsData = b;                 }

    void                SetProcessID(Int_t n)             { fDecay = n;                  }
    void                SetTupleName(const char* c)       { fTupleName = c;              }
    void                SetFillNtupleType(Int_t d)        { fFillNtupleType= d;          }

  protected:
    void                Process();
    void                SlaveBegin();
    void                SlaveTerminate();
    // Private auxiliary methods...
    // Names of the input Collections
    TString             fMetName;
    TString             fPhotonsName;
    TString             fElectronsName;
    TString             fMuonsName;
    TString             fTausName;
    TString             fJetsName;
    TString             fLeptonsName;
    TString             fVertexName;

    TString             fSuperClustersName;
    TString             fTracksName;
    TString             fPVName;
    TString             fPileUpDenName;    
    TString             fPileUpName;
    TString             fBeamspotName;
    TString             fMCEvInfoName;
    TString             fMCPartName;

    // is it Data or MC?
    Bool_t              fIsData;
    
    // there are not some PV pre-selection?
    Bool_t              fMetFromBranch;
    Bool_t              fPhotonsFromBranch;
    Bool_t              fElectronsFromBranch;
    Bool_t              fMuonsFromBranch;
    Bool_t              fTausFromBranch;
    Bool_t              fJetsFromBranch;
    Bool_t              fPVFromBranch;
    Bool_t              fQGTaggerCHS;

    QGTagger           *qgTagger;

    const PFMetCol                *fRawMet;
    const MetCol                  *fMet;
    const PhotonCol               *fPhotons;
    const ElectronCol             *fElectrons;
    const MuonCol                 *fMuons;
    const PFTauCol                *fPFTaus;
    const JetCol                  *fJets;
    const EvtSelData              *fEvtSelData;

    const TrackCol                *fTracks;
    const VertexCol               *fPV;
    const BeamSpotCol             *fBeamspot;
    const MCEventInfo             *fMCEventInfo;
    const PileupInfoCol           *fPileUp;    
    const PileupEnergyDensityCol  *fPileUpDen;
    const SuperClusterCol         *fSuperClusters; 
    const MCParticleCol           *fParticles;	        

    // --------------------------------
    Int_t                          fDecay;
    TFile	                      *fOutputFile;
    TString	                       fTupleName;
    Int_t                          fFillNtupleType;
    MitGPTree                      fMitGPTree;

    Int_t                          fNEventsSelected;

    ClassDef(MonoJetTreeWriter, 1) // Photon identification module
      };
}
#endif

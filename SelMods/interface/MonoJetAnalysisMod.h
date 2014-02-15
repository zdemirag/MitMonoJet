//--------------------------------------------------------------------------------------------------
// MonoJetAnalysisMod
//
// An analysis module for selecting gamma+MET events and produces some basic distributions.
//
//
// Authors: LDM, TJW
//
//--------------------------------------------------------------------------------------------------
#ifndef MITMonoJet_SELMODS_MonoJetANALYSISMOD_H
#define MITMonoJet_SELMODS_MonoJetANALYSISMOD_H

#include "MitAna/TreeMod/interface/BaseMod.h" 
#include "MitAna/DataTree/interface/MetCol.h"
#include "MitAna/DataTree/interface/ElectronCol.h"
#include "MitAna/DataTree/interface/MuonCol.h"
#include "MitAna/DataTree/interface/PFTauCol.h"
#include "MitAna/DataTree/interface/JetCol.h"
#include "MitAna/DataTree/interface/CaloJetCol.h"
#include "MitAna/DataTree/interface/PFJetCol.h"
#include "MitAna/DataTree/interface/SuperClusterCol.h"
#include "MitAna/DataTree/interface/DecayParticleCol.h"
#include "MitAna/DataTree/interface/PFCandidateCol.h"
#include "MitAna/DataTree/interface/VertexCol.h"
#include "MitPhysics/Utils/interface/MuonTools.h"
#include "MitAna/DataTree/interface/PFMetCol.h"
#include "MitAna/DataTree/interface/PhotonCol.h"

class TH1D;
class TH2D;

namespace mithep 
{
  class MonoJetAnalysisMod : public BaseMod
  {
  public:
    MonoJetAnalysisMod(const char *name  = "MonoJetAnalysisMod", 
		       const char *title = "MoneJet Slection Module");
    ~MonoJetAnalysisMod() {}
    
    // setting all the input Names
    void                  SetInputMetName          (const char *n) { fMetBranchName= n;           }
    void 		  SetJetsName              (const char *n) { fJetsName = n;               } 
    void                  SetElectronsName         (const char *n) { fElectronsName = n;          }
    void                  SetMuonsName             (const char *n) { fMuonsName = n;              }
    void                  SetTausName              (const char *n) { fTausName = n;               }
    void                  SetLeptonsName           (const char *n) { fLeptonsName = n;            }
    void                  SetPFCandidatesName      (const char *n) { fPFCandidatesName = n;       }
    void                  SetPFCandidatesFromBranch(bool b)        { fPFCandidatesFromBranch = b; }

    // decide whether to read from branch
    void                  SetMetFromBranch         (Bool_t b)      { fMetFromBranch = b;          }
    void                  SetJetsFromBranch        (Bool_t b)      { fJetsFromBranch = b;         }
    void                  SetElectronsFromBranch   (Bool_t b)      { fElectronsFromBranch = b;    }
    void                  SetMuonsFromBranch       (Bool_t b)      { fMuonsFromBranch = b;        }
    void                  SetTausFromBranch        (Bool_t b)      { fTausFromBranch = b;         }
    
  protected:
    // Standard module methods
    void                  Begin();
    void                  Process();
    void                  SlaveBegin();
    void                  SlaveTerminate();
    void                  Terminate();      
    // Setting cut values
    void                  SetMinNumLeptons         (Int_t n)       { fMinNumLeptons = n;          }
    void                  SetMinNumTaus            (Int_t n)       { fMinNumTaus = n;             }
    void 	          SetMinNumJets            (Int_t n)       { fMinNumJets = n;             }
    void 	          SetMinJetEt              (Double_t x)    { fMinJetEt = x;               }
    void	          SetMaxJetEta             (Double_t x)    { fMaxJetEta = x;              }
    void 	          SetMinMetEt              (Double_t x)    { fMinMetEt = x;               }
    void                  SetMinChargedHadronFrac  (Double_t x)    { fMinChargedHadronFrac = x;   }
    void                  SetMaxNeutralHadronFrac  (Double_t x)    { fMaxNeutralHadronFrac = x;   }
    void                  SetMaxNeutralEmFrac      (Double_t x)    { fMaxNeutralEmFrac = x;       }

    // names of the collections
    TString               fMetBranchName;
    TString		  fJetsName;
    TString               fElectronsName;
    TString               fMuonsName;
    TString               fTausName;
    TString               fLeptonsName;
    TString               fPFCandidatesName;
    // logical whether to read from branch
    Bool_t                fMetFromBranch;
    Bool_t                fJetsFromBranch;
    Bool_t                fPFCandidatesFromBranch;
    Bool_t                fElectronsFromBranch;
    Bool_t                fMuonsFromBranch;
    Bool_t                fTausFromBranch;
    // hooks to the collections
    const MetCol         *fMet;
    const JetCol	 *fJets; 
    const ElectronCol    *fElectrons;
    const MuonCol        *fMuons;
    const PFTauCol       *fPFTaus;
    const PFCandidateCol *fPFCandidates;

    // Cuts
    UInt_t                fMinNumLeptons;
    UInt_t                fMinNumTaus;
    UInt_t                fMinNumJets;
    Double_t              fMinJetEt;
    Double_t              fMaxJetEta;
    Double_t              fMinMetEt;
    Double_t              fMinChargedHadronFrac;
    Double_t              fMaxNeutralHadronFrac;
    Double_t              fMaxNeutralEmFrac;
    
    // Counters
    Int_t                 fNEventsSelected;

    // Histograms
    TH1D                 *fMonoJetSelection;         // cut flow monitoring
    TH1D                 *fPhotonEt;                 // photon transverse energy spectrum
    TH1D                 *fMetEt;                    // met spectrum
    TH1D		 *fJetEt;       	     // jet Et spectrum
    TH1D		 *fJetEta;                   // jet eta
    
    ClassDef(MonoJetAnalysisMod,1) // MonJet Selection Module
  };
}
#endif

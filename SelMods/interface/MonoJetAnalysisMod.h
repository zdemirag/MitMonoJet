//--------------------------------------------------------------------------------------------------
// MonoJetAnalysisMod
//
// An analysis module for selecting gamma+MET events and produces some basic distributions.
//
//
// Authors: LDM, TJW, YI
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
#include "MitAna/DataTree/interface/TriggerMask.h"

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
    void                  SetCategoriesName        (const char *n) { fCategoriesName = n; }
    void                  SetInputMetName          (const char *n) { fMetBranchName= n;           }
    void                  SetJetsName              (const char *n) { fJetsName = n;               } 
    void                  SetElectronsName         (const char *n) { fElectronsName = n;          }
    void                  SetMuonsName             (const char *n) { fMuonsName = n;              }
    void                  SetTausName              (const char *n) { fTausName = n;               }
    void                  SetLeptonsName           (const char *n) { fLeptonsName = n;            }
    void                  SetPFCandidatesName      (const char *n) { fPFCandidatesName = n;       }
    void                  SetPFCandidatesFromBranch(bool b)        { fPFCandidatesFromBranch = b; }

    char const*           GetCategoriesName() const { return fCategoriesName.Data(); }

    // decide whether to read from branch
    void                  SetMetFromBranch         (Bool_t b)      { fMetFromBranch = b;          }
    void                  SetJetsFromBranch        (Bool_t b)      { fJetsFromBranch = b;         }
    void                  SetElectronsFromBranch   (Bool_t b)      { fElectronsFromBranch = b;    }
    void                  SetMuonsFromBranch       (Bool_t b)      { fMuonsFromBranch = b;        }
    void                  SetTausFromBranch        (Bool_t b)      { fTausFromBranch = b;         }

    // Setting cut values
    void                  SetMinNumLeptons         (UInt_t c, Int_t n)       { fMinNumLeptons[c] = n;          }
    void                  SetMaxNumLeptons         (UInt_t c, Int_t n)       { fMaxNumLeptons[c] = n;          }
    void                  SetMinNumTaus            (UInt_t c, Int_t n)       { fMinNumTaus[c] = n;             }
    void                  SetMaxNumTaus            (UInt_t c, Int_t n)       { fMaxNumTaus[c] = n;             }
    void                  SetMinNumJets            (UInt_t c, Int_t n)       { fMinNumJets[c] = n;             }
    void                  SetMaxNumJets            (UInt_t c, Int_t n)       { fMaxNumJets[c] = n;             }
    void                  SetMinNumGenNeutrinos    (UInt_t c, Int_t n)       { fMinNumGenNeutrinos[c] = n;     }
    void                  SetMinJetEt              (UInt_t c, Double_t x)    { fMinJetEt[c] = x;               }
    void                  SetMaxJetEta             (UInt_t c, Double_t x)    { fMaxJetEta[c] = x;              }
    void                  SetMinMetEt              (UInt_t c, Double_t x)    { fMinMetEt[c] = x;               }
    void                  SetMinEmulMetEt          (UInt_t c, Double_t x)    { fMinEmulMetEt[c] = x;           }
    void                  SetMinChargedHadronFrac  (UInt_t c, Double_t x)    { fMinChargedHadronFrac[c] = x;   }
    void                  SetMaxNeutralHadronFrac  (UInt_t c, Double_t x)    { fMaxNeutralHadronFrac[c] = x;   }
    void                  SetMaxNeutralEmFrac      (UInt_t c, Double_t x)    { fMaxNeutralEmFrac[c] = x;       }
    
  protected:
    // Standard module methods
    void                  Begin();
    void                  Process();
    void                  SlaveBegin();
    void                  SlaveTerminate();
    void                  Terminate();      

    // names of the collections
    TString               fCategoriesName;
    TString               fMetBranchName;
    TString               fJetsName;
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
    const JetCol         *fJets; 
    const ElectronCol    *fElectrons;
    const MuonCol        *fMuons;
    const PFTauCol       *fPFTaus;
    const PFCandidateCol *fPFCandidates;

    static UInt_t const   nCat = 1 * 8; // must be a multiple of 8

    // Cuts
    // The module can define up to nCat set of cuts, each corresponing
    // to a skim category.
    // Pretty bad implementation as the output file does not document
    // what cuts resulted in which cateogry, but this is just an example..
    UInt_t                fMinNumLeptons[nCat];
    UInt_t                fMaxNumLeptons[nCat];
    UInt_t                fMinNumTaus[nCat];
    UInt_t                fMaxNumTaus[nCat];
    UInt_t                fMinNumJets[nCat];
    UInt_t                fMaxNumJets[nCat];
    UInt_t                fMinNumGenNeutrinos[nCat];
    Double_t              fMinJetEt[nCat];
    Double_t              fMaxJetEta[nCat];
    Double_t              fMinMetEt[nCat];
    Double_t              fMinEmulMetEt[nCat];
    Double_t              fMinChargedHadronFrac[nCat];
    Double_t              fMaxNeutralHadronFrac[nCat];
    Double_t              fMaxNeutralEmFrac[nCat];

    // Currently we do not have a generic named trigger bits object
    // Must encapsulate the TriggerMask in ObjArray (with a single element)
    mithep::Array<mithep::TriggerMask>*    fCategories;
    
    // Counters
    Int_t                 fNEventsSelected;

    // Histograms
    TH1D                 *fPhotonEt[nCat];                 // photon transverse energy spectrum
    TH1D                 *fMetEt[nCat];                    // met spectrum
    TH1D                 *fJetEt[nCat];                // jet Et spectrum
    TH1D                 *fJetEta[nCat];                   // jet eta
    
    ClassDef(MonoJetAnalysisMod,1) // MonJet Selection Module
  };
}
#endif

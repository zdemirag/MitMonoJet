//--------------------------------------------------------------------------------------------------
// BoostedVAnalysisMod
//
// An analysis module for selecting jet+MET events.
//
//
// Authors: LDM
//
//--------------------------------------------------------------------------------------------------
#ifndef MITMONOJET_SELMODS_BOOSTEDVANALYSISMOD_H
#define MITMONOJET_SELMODS_BOOSTEDVANALYSISMOD_H

#include "MitAna/TreeMod/interface/BaseMod.h" 
#include "MitAna/DataTree/interface/MetCol.h"
#include "MitAna/DataTree/interface/ElectronCol.h"
#include "MitAna/DataTree/interface/MuonCol.h"
#include "MitAna/DataTree/interface/JetCol.h"
#include "MitPhysics/Utils/interface/MuonTools.h"
#include "MitAna/DataTree/interface/PFMetCol.h"

class TH1D;
class TH2D;

namespace mithep 
{
  class BoostedVAnalysisMod : public BaseMod
  {
  public:
    BoostedVAnalysisMod(const char *name  = "BoostedVAnalysisMod", 
                        const char *title = "Boosted V selection Module");
    ~BoostedVAnalysisMod() {}
    
    // setting all the input Names
    void                  SetInputMetName          (const char *n) { fMetBranchName= n;           }
    void                  SetJetsName              (const char *n) { fJetsName = n;               } 
    void                  SetElectronsName         (const char *n) { fElectronsName = n;          }
    void                  SetMuonsName             (const char *n) { fMuonsName = n;              }
    void                  SetLeptonsName           (const char *n) { fLeptonsName = n;            }

    // decide whether to read from branch
    void                  SetMetFromBranch         (Bool_t b)      { fMetFromBranch = b;          }
    void                  SetJetsFromBranch        (Bool_t b)      { fJetsFromBranch = b;         }
    void                  SetElectronsFromBranch   (Bool_t b)      { fElectronsFromBranch = b;    }
    void                  SetMuonsFromBranch       (Bool_t b)      { fMuonsFromBranch = b;        }

    // decide whether to apply given preselections or not
    void                  ApplyTopPresel           (Bool_t b)      { fApplyTopPresel = b;         }
    void                  ApplyWlepPresel          (Bool_t b)      { fApplyWlepPresel = b;        }
    void                  ApplyZlepPresel          (Bool_t b)      { fApplyZlepPresel = b;        }
    void                  ApplyMetPresel           (Bool_t b)      { fApplyMetPresel = b;         }

    // Setting cut values
    void                  SetMinTagJetPt           (Double_t x)    { fMinTagJetPt = x;            }
    void                  SetMinMet                (Double_t x)    { fMinMet = x;                 }
    
  protected:
    // Standard module methods
    void                  Begin();
    void                  Process();
    void                  SlaveBegin();
    void                  SlaveTerminate();
    void                  Terminate();      
    // Auxiliary methods
    void                  CorrectMet(const float met, const float metPhi,
                                     const Particle *l1, const Particle *l2,
                                     float &newMet, float &newMetPhi);
    // names of the collections
    TString               fMetBranchName;
    TString               fJetsName;
    TString               fElectronsName;
    TString               fMuonsName;
    TString               fLeptonsName;
    // logical whether to read from branch
    Bool_t                fMetFromBranch;
    Bool_t                fJetsFromBranch;
    Bool_t                fElectronsFromBranch;
    Bool_t                fMuonsFromBranch;
    // logical whether to apply given preselections or not
    Bool_t                fApplyTopPresel; 
    Bool_t                fApplyWlepPresel;
    Bool_t                fApplyZlepPresel;
    Bool_t                fApplyMetPresel; 
    // hooks to the collections
    const PFMetCol       *fMet;
    const JetCol         *fJets; 
    const ElectronCol    *fElectrons;
    const MuonCol        *fMuons;
    // Cuts
    Double_t              fMinTagJetPt;
    Double_t              fMinMet;
    // Counters
    Long64_t              fAll;
    Long64_t              fPass;
        
    ClassDef(BoostedVAnalysisMod,1) // MonJet Selection Module
  };
}
#endif

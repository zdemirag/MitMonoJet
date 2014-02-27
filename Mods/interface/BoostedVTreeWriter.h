//--------------------------------------------------------------------------------------------------
// BoostedVTreeWriter
//
// This module writes a tree about jets corresponding to Vector bosons. See documentation here:
//
// Authors: C.Paus
//--------------------------------------------------------------------------------------------------
#ifndef MITMONOJET_MODS_BoostedVTreeWriter_H
#define MITMONOJET_MODS_BoostedVTreeWriter_H

#include <TH1.h>

#include "fastjet/JetDefinition.hh"
#include "fastjet/GhostedAreaSpec.hh"
#include "fastjet/AreaDefinition.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/tools/Pruner.hh"

#include "MitAna/TreeMod/interface/BaseMod.h"
#include "MitAna/DataTree/interface/TriggerObjectCol.h"
#include "MitAna/DataTree/interface/JetCol.h"
#include "MitAna/DataTree/interface/PFJetCol.h"
#include "MitAna/DataTree/interface/PFCandidateCol.h"
#include "MitMonoJet/Core/MitGPBoostedVTree.h"


namespace mithep
{
  class BoostedVTreeWriter : public BaseMod
  {
  public:
    BoostedVTreeWriter(const char *name  = "BoostedVTreeWriter",
		       const char *title = "Vector Boson Tagging module");

    void                          SetTriggerObjsName(const char *n) { fTriggerObjsName = n;   }
    void                          SetHistNPtBins(Int_t n)           { fHistNPtBins = n;       }
    void                          SetHistNEtaBins(Int_t n)          { fHistNEtaBins = n;      }
    void                          SetHistMinPt(Double_t b)          { fHistMinPt = b;         }
    void                          SetHistMaxPt(Double_t b)          { fHistMaxPt = b;         }
    void                          SetHistMinEta(Double_t b)         { fHistMinEta = b;        }
    void                          SetHistMaxEta(Double_t b)         { fHistMaxEta = b;        }
    void                          SetHistTau1Bins(Int_t n)          { fHistTau1Bins = n;      }
    void                          SetHistTau2Bins(Int_t n)          { fHistTau2Bins = n;      }
    void                          SetHistTau3Bins(Int_t n)          { fHistTau3Bins = n;      }
    void                          SetHistT2ovrT1Bins(Double_t b)    { fHistT2ovrT1Bins = b;   }
    void                          SetHistT3ovrT2Bins(Double_t b)    { fHistT3ovrT2Bins = b;   }
    void                          SetHistMinTau1(Double_t b)        { fHistMinTau1 = b;       }
    void                          SetHistMinTau2(Double_t b)        { fHistMinTau2 = b;       }
    void                          SetHistMinTau3(Double_t b)        { fHistMinTau3 = b;       }
    void                          SetHistMinT2ovrT1(Double_t b)     { fHistMinT2ovrT1 = b;    }
    void                          SetHistMinT3ovrT2(Double_t b)     { fHistMinT3ovrT2 = b;    }
    void                          SetHistMaxTau1(Double_t b)        { fHistMaxTau1 = b;       }
    void                          SetHistMaxTau2(Double_t b)        { fHistMaxTau2 = b;       }
    void                          SetHistMaxTau3(Double_t b)        { fHistMaxTau3 = b;       }
    void                          SetHistMaxT2ovrT1(Double_t b)     { fHistMaxT2ovrT1 = b;    }
    void                          SetHistMaxT3ovrT2(Double_t b)     { fHistMaxT3ovrT2 = b;    }

  private:
    float                         GetTau(fastjet::PseudoJet &iJet,int iN, float iKappa);

    TString                       fTriggerObjsName;       //(i) name of trigger objects
    const TriggerObjectCol       *fTrigObjs;              //trigger objects coll handle
    TString                       fJetsName;              //(i) name of jets used to make trigger
    const JetCol                 *fJets;                  //jets used to make the trigger
    TString                       fPFCandidatesName;      //(i) name of PF candidates coll
    const PFCandidateCol         *fPFCandidates;          //particle flow candidates coll handle

    // Objects from fastjet we want to use
    double                        fConeSize;
    fastjet::Pruner              *fPruner;
    fastjet::JetDefinition       *fCAJetDef;
    fastjet::GhostedAreaSpec     *fActiveArea;
    fastjet::AreaDefinition      *fAreaDefinition;

  protected:
    void                          Process();
    void                          SlaveBegin();
    void                          Terminate();
    // Output histograms
    TH1D                         *fPFCandidatesPt;        //particle flow pt  distribution
    TH1D                         *fPFCandidatesEta;       //particle flow eta distribution
    TH1D                         *fCAJetPt;               //fastjet output pt distribution
    TH1D                         *fCAJetEta;              //fastjet output eta distribution
    TH1D                         *fCATau1;                //plot of Tau 1
    TH1D                         *fCATau2;                //plot of Tau 2
    TH1D                         *fCATau3;                //plot of Tau 3
    TH1D                         *fCAT2ovrT1;             //plot of Tau 2 over Tau 1
    TH1D                         *fCAT3ovrT2;             //plot of Tau 3 over Tau 2

    Int_t                         fHistNPtBins;
    Int_t                         fHistNEtaBins;
    Double_t                      fHistMinPt;
    Double_t                      fHistMaxPt;
    Double_t                      fHistMinEta;
    Double_t                      fHistMaxEta;
    Int_t                         fHistTau1Bins;
    Int_t                         fHistTau2Bins;
    Int_t                         fHistTau3Bins;
    Int_t                         fHistT2ovrT1Bins;
    Int_t                         fHistT3ovrT2Bins;
    Double_t                      fHistMinTau1;
    Double_t                      fHistMinTau2;
    Double_t                      fHistMinTau3;
    Double_t                      fHistMinT2ovrT1;
    Double_t                      fHistMinT3ovrT2;
    Double_t                      fHistMaxTau1;
    Double_t                      fHistMaxTau2;
    Double_t                      fHistMaxTau3;
    Double_t                      fHistMaxT2ovrT1;
    Double_t                      fHistMaxT3ovrT2;

    // Output tree
    MitGPBoostedVTree             fMitGPTree;

    ClassDef(BoostedVTreeWriter, 0) // Boosted Vector boson tree writer
  };
}
#endif

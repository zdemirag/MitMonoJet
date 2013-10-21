//--------------------------------------------------------------------------------------------------
// $Id $
//
// MonoJetAnalysisMod
//
// A Module for Selecting gamma+MET events
// and produces some distributions.
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
    MonoJetAnalysisMod(const char *name="MonoJetAnalysisMod", 
		 const char *title="Example analysis module with all branches");
      ~MonoJetAnalysisMod() {}

      // setting all the input Names
      void                SetInputMetName    (const char *n){ fMetBranchName= n; }
      void 		        SetJetsName        (const char *n){ fJetsName = n;     } 
      void                SetElectronsName   (const char *n){ fElectronsName = n;}
      void                SetMuonsName       (const char *n){ fMuonsName = n;    }
      void                SetTausName        (const char *n){ fTausName = n;     }
      void                SetLeptonsName     (const char *n){ fLeptonsName = n;  }

      void                SetMetFromBranch(Bool_t b)        { fMetFromBranch = b;      }
      void                SetJetsFromBranch(Bool_t b)       { fJetsFromBranch = b;     }
      void                SetElectronsFromBranch(Bool_t b)  { fElectronsFromBranch = b;}
      void                SetMuonsFromBranch(Bool_t b)      { fMuonsFromBranch = b;    }
      void                SetTausFromBranch(Bool_t b)       { fTausFromBranch = b;     }

    protected:
      TString                  fMetBranchName;           //name of input met branch
      TString		       fJetsName;                //name of input jet branch (added by TJ)
      TString                  fElectronsName;
      TString                  fMuonsName;
      TString                  fTausName;
      TString                  fLeptonsName;

      Bool_t                   fMetFromBranch;           //met is loaded from a branch
      Bool_t                   fJetsFromBranch;          //jet are loaded from a branch
      Bool_t                   fElectronsFromBranch;
      Bool_t                   fMuonsFromBranch;
      Bool_t                   fTausFromBranch;

      TH1D                    *fMonoJetSelection;        //histogram for cut flow monitoring
      TH1D                    *fPhotonEt;                //histogram of photon transverse energy spectrum
      TH1D                    *fMetEt;                   //histogram of met spectrum
      TH1D		            *fJetEt;       	         //histogram of jet spectrum (added by TJ)
      TH1D		            *fJetEta;			   //histogram of jet eta (added by TJ; for testing purposes)

      const MetCol            *fMet;
      const JetCol	      *fJets; 
      const ElectronCol       *fElectrons;
      const MuonCol           *fMuons;
      const PFTauCol          *fPFTaus;
      
      void     SetMinNumLeptons(Int_t n)  { fMinNumLeptons = n; }
      void     SetMinNumTaus(Int_t n)     { fMinNumTaus = n;    }
      void 	   SetMinNumJets(Int_t n)     { fMinNumJets = n;    }
      void 	   SetMinJetEt(Double_t x)    { fMinJetEt = x;      }
      void	   SetMaxJetEta(Double_t x)   { fMaxJetEta = x;     }
      void 	   SetMinMetEt(Double_t x)    { fMinMetEt = x;      }
      void     SetMinChargedHadronFrac(Double_t x) { fMinChargedHadronFrac = x; }
      void     SetMaxNeutralHadronFrac(Double_t x) { fMaxNeutralHadronFrac = x; }
      void     SetMaxNeutralEmFrac(Double_t x)     { fMaxNeutralEmFrac = x;     }

      void         Begin();
      void         Process();
      void         SlaveBegin();
      void         SlaveTerminate();
      void         Terminate();      

      unsigned int fMinNumLeptons;
      unsigned int fMinNumTaus;
      unsigned int fMinNumJets;
      Double_t     fMinJetEt;
      Double_t     fMaxJetEta;
      Double_t     fMinMetEt;
      Double_t     fMinChargedHadronFrac;
      Double_t     fMaxNeutralHadronFrac;
      Double_t     fMaxNeutralEmFrac;

      Int_t        fNEventsSelected; //selected events

      ClassDef(MonoJetAnalysisMod,1) // TAM example analysis module
  };
}
#endif

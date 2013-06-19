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
      void 		          SetJetsName        (const char *n){ fJetsName = n; } //added by TJ
      void                SetMetFromBranch(Bool_t b)    { fMetFromBranch = b; }
      void                SetJetsFromBranch(Bool_t b)    { fJetsFromBranch = b; }
    protected:
      TString                  fMetBranchName;           //name of input met branch
      TString		           fJetsName;              	 //name of input jet branch (added by TJ)
      Bool_t                   fMetFromBranch;           //met is loaded from a branch
      Bool_t                   fJetsFromBranch;          //jet are loaded from a branch

      TH1D                    *fMonoJetSelection;        //histogram for cut flow monitoring
      TH1D                    *fPhotonEt;                //histogram of photon transverse energy spectrum
      TH1D                    *fMetEt;                   //histogram of met spectrum
      TH1D		              *fJetEt;		          	 //histogram of jet spectrum (added by TJ)
      TH1D		              *fJetEta;			         //histogram of jet eta (added by TJ; for testing purposes)

      const PFMetCol          *fMet;
      const JetCol		      *fJets; //added by TJ

      void 	   SetMinNumJets(Int_t n)  { fMinNumJets = n; }
      void 	   SetMinJetEt(Double_t x) { fMinJetEt = x; }
      void	   SetMaxJetEta(Double_t x){ fMaxJetEta = x; }
      void 	   SetMinMetEt(Double_t x) { fMinMetEt = x; }

      void         Begin();
      void         Process();
      void         SlaveBegin();
      void         SlaveTerminate();
      void         Terminate();      

      unsigned int fMinNumJets;
      Double_t fMinJetEt;
      Double_t fMaxJetEta;
      Double_t fMinMetEt;

      Int_t                    fNEventsSelected;         //selected events

      ClassDef(MonoJetAnalysisMod,1) // TAM example analysis module
  };
}
#endif

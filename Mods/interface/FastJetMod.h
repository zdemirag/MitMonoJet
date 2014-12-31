//--------------------------------------------------------------------------------------------------
// $Id: FastJetMod.h,v 1.9 2011/03/01 17:27:22 mzanetti Exp $
//
// FastJetMod 
//
// This module process a collection of input PFCandidates, cluster them
// using fastjet and then spits out Bambu PFJets collection(s)
//
// Authors: L.DiMatteo
//--------------------------------------------------------------------------------------------------

#ifndef MITMONOJET_MODS_FastJetMod_H
#define MITMONOJET_MODS_FastJetMod_H

#include "fastjet/PseudoJet.hh"
#include "fastjet/JetDefinition.hh"
#include "fastjet/GhostedAreaSpec.hh"
#include "fastjet/AreaDefinition.hh"
#include "fastjet/ClusterSequenceArea.hh"

#include "MitAna/TreeMod/interface/BaseMod.h"
#include "MitAna/DataTree/interface/JetCol.h"
#include "MitAna/DataTree/interface/PFJetCol.h"
#include "MitAna/DataTree/interface/PFCandidateCol.h"

namespace mithep
{
  class FastJetMod : public BaseMod
  {
    public:
      FastJetMod(const char *name = "FastJetMod",
                   const char *title = "FastJet module");
      ~FastJetMod();

      const char *GetOutputJetsName()     const { return fOutputJetsName;     }

      void GetMatchBtag(Bool_t b)               { fGetMatchBtag = b;          }
      void UseBambuJets(Bool_t b)               { fUseBambuJets = b;          }

      void SetBtaggedJetsName(const char *n)    { fBtaggedJetsName = n;       }
      void SetBtaggedJetsFromBranch(Bool_t b)   { fBtaggedJetsFromBranch = b; }
                                                                                                                                                  
      void SetJetsName(const char *n)           { fJetsName = n;              }
      void SetJetsFromBranch(Bool_t b)          { fJetsFromBranch = b;        }
                                                                              
      void SetPfCandidatesName(const char *n)   { fPfCandidatesName = n;      }
      void SetPfCandidatesFromBranch(Bool_t b)  { fPfCandidatesFromBranch = b;}

      void SetOutputJetsName(const char *n)     { fOutputJetsName = n;        }
                                                                    
      void SetConeSize(double d)                { fJetConeSize = d;           }

      void SetParticleMinPt(double d)           { fParticleMinPt = d;         }
      void SetJetMinPt(double d)                { fJetMinPt = d;              }
 
    protected:
      void Process();
      void SlaveBegin();
      void SlaveTerminate();
      
      // PFJet filler helper
      void FillPFJet (PFJet *pPFJet, fastjet::PseudoJet &fjJet);    
        
    private:
      Bool_t fGetMatchBtag;                //=true if b-tag obtained by match with standard jets (AK5)
      Bool_t fUseBambuJets;                //=true if input small jets already present in bambu
      Bool_t fUseBambuFatJets;             //=true if input large jets already present in bambu

      TString fBtaggedJetsName;            //(i) name of input btagged jets
      Bool_t fBtaggedJetsFromBranch;       //are input btagged jets from Branch?
      const PFJetCol *fBtaggedJets;        //input btagged jets
      
      TString fJetsName;                   //(i) name of input jets
      Bool_t fJetsFromBranch;              //are input jets from Branch?
      const PFJetCol *fJets;               //input jets

      TString fPfCandidatesName;           //(i) name of PF candidates coll
      Bool_t fPfCandidatesFromBranch;      //are PF candidates from Branch?
      const PFCandidateCol *fPfCandidates; //particle flow candidates coll handle
 
      TString fOutputJetsName;             //name of output jets collection
      JetOArr *fOutputJets;                //output jets collection
      
      // Objects from fastjet we want to use
      double fJetConeSize;                 //fastjet clustering radius
      double fFatJetConeSize;              //fastjet fat clustering radius
      fastjet::JetDefinition *fJetDef  ;   //fastjet clustering definition
      fastjet::GhostedAreaSpec *fActiveArea;
      fastjet::AreaDefinition *fAreaDefinition;

      double fParticleMinPt;               //low cutoff to avoid including unphysical PFCands
      double fJetMinPt;                    //low cutoff to carrying on soft jets
      
      ClassDef(FastJetMod, 0)              //FastJet bambu producer      
  };
}
#endif

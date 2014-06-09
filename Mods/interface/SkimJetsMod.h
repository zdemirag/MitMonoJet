//--------------------------------------------------------------------------------------------------
// SkimJetsMod
//
// This module is derived from the templated SkimJetsMod contained in MitAna/PhysicsMod
// Its developement is due to the fact that cleaned jets are saved as Jets while 
// they are in fact created as PFJets. A dynamic conversion is performed here in 
// order to store the correct classtype in the output collection
//
// Authors: L.DiMatteo, C.Paus
//--------------------------------------------------------------------------------------------------
#ifndef MITMONOJET_MODS_SKIMJETSMOD_H
#define MITMONOJET_MODS_SKIMJETSMOD_H

#include "MitAna/DataCont/interface/Array.h"
#include "MitAna/DataCont/interface/ObjArray.h"
#include "MitAna/DataTree/interface/JetCol.h"
#include "MitAna/DataTree/interface/PFJetCol.h"
#include "MitAna/TreeMod/interface/BaseMod.h" 

namespace mithep 
{
  class SkimJetsMod : public BaseMod 
  {
    public:
      SkimJetsMod(const char *name  = "SkimJetsMod", 
	      const char *title = "Jet skim module");
      ~SkimJetsMod();

      const char              *GetBranchName()              const { return fBranchName; }
      void                     SetBranchName(const char *n)       { fBranchName = n;    }
      void                     SetColFromBranch(bool b)           { fColFromBranch = b; }
      void                     SetColMarkFilter(bool b)           { fColMarkFilter = b; }
      void                     SetPublishArray(bool b)            { fPublishArray = b;  }

    protected:
      TString                  fBranchName;    //name of collection
      bool                     fColFromBranch; //input collection is branch?
      bool                     fColMarkFilter; //input collection filter based on mark?
      bool                     fPublishArray;  //output collection is Array?
  
      const Collection<Jet>    *fCol;           //!pointer to collection (in - branch) 
      ObjArray<PFJet>          *fColSkm;        //!pointer to collection (skm - published object)
      Array<PFJet>             *fArrSkm;        //!pointer to array (skm - published object)

      void                     Process();
      void                     SlaveBegin();
      void                     SlaveTerminate();

      ClassDef(SkimJetsMod, 0) // Skim jets module
  };
}
#endif

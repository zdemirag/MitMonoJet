//--------------------------------------------------------------------------------------------------
// MetAnalysisMod 
//
// A minimal module to plot the missing Et
//
// Authors: Y.Iiyama
//--------------------------------------------------------------------------------------------------

#ifndef MITMONOJET_MODS_MetAnalysisMod_H
#define MITMONOJET_MODS_MetAnalysisMod_H

#include "MitAna/TreeMod/interface/BaseMod.h"
#include "MitAna/DataTree/interface/MetFwd.h"
#include "MitAna/DataTree/interface/TriggerMask.h"
#include "MitAna/DataCont/interface/Array.h"

#include "TH1D.h"

namespace mithep
{
  class MetAnalysisMod : public BaseMod
  {
    public:
      MetAnalysisMod(char const* name = "MetAnalysisMod",
                     char const* title = "Met Analysis module");
      ~MetAnalysisMod();

      void SetMetName(char const* n) { fMetColName = n; }
      void SetMetColFromBranch(bool b) { fMetColFromBranch = b; }
      void SetMonoJetCategoriesName(char const* n) { fMonoJetCategoriesName = n; }
      void SetMonoJetCategory(int cat) { fMonoJetCategory = cat; }
 
    protected:
      void Process();
      void SlaveBegin();
      void SlaveTerminate();
      
    private:
      TString fMetColName;
      bool fMetColFromBranch;
      MetCol const* fMetCol;
      int fMonoJetCategory;
      TString fMonoJetCategoriesName;
      mithep::Array<mithep::TriggerMask> const* fMonoJetCategories;
      
      TH1D* hMet;
      
      ClassDef(MetAnalysisMod, 1) //Met analysis module
  };
}
#endif

#include "MitMonoJet/Mods/interface/MetAnalysisMod.h"
#include "MitAna/DataTree/interface/Met.h"

ClassImp(mithep::MetAnalysisMod)

namespace mithep {

  MetAnalysisMod::MetAnalysisMod(const char *name/* = "MetAnalysisMod"*/,
                                 const char *title/* = "Met Analysis module"*/) :
    BaseMod(name, title),
    fMetColName("PFMet"),
    fMetCol(0),
    fMonoJetCategory(-1),
    fMonoJetCategoriesName("MonoJetEventCategories"),
    fMonoJetCategories(0),
    hMet(0)
  {
  }

  MetAnalysisMod::~MetAnalysisMod()
  {
  }

  void
  MetAnalysisMod::Process()
  {
    // if a mono-jet skim category is set, only use events in that category
    if(fMonoJetCategory >= 0){
      if(!LoadEventObject(fMonoJetCategoriesName, fMonoJetCategories, true))
        return;

      if(!fMonoJetCategories->At(0)->TestBit(fMonoJetCategory))
        return;
    }

    if(!LoadEventObject(fMetColName, fMetCol, true))
      return;

    mithep::Met const& met(*fMetCol->At(0));

    hMet->Fill(met.Pt());
  }

  void
  MetAnalysisMod::SlaveBegin()
  {
    AddTH1(hMet, "hMet", "Missing E_{T};E_{T}^{miss} (GeV);events / 50 GeV", 20, 0., 1000.);

    ReqEventObject(fMetColName, fMetCol, fMetColFromBranch);
    if(fMonoJetCategory >= 0)
      ReqEventObject(fMonoJetCategoriesName, fMonoJetCategories, true);
  }

  void
  MetAnalysisMod::SlaveTerminate()
  {
  }

}

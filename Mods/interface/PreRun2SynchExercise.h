#ifndef MITMONOJET_MODS_PRERUN2SYNCHEXERCISE_H
#define MITMONOJET_MODS_PRERUN2SYNCHEXERCISE_H

#include "MitAna/TreeMod/interface/BaseMod.h"

#include "TString.h"
#include "TTree.h"

namespace mithep {

  class PreRun2SynchExercise : public BaseMod {
  public:
    PreRun2SynchExercise(char const* name = "PreRun2SynchExercise", char const* title = "Pre-Run2 synch") : BaseMod(name, title) {}
    ~PreRun2SynchExercise() {}

  protected:
    void Process() override;
    void SlaveBegin() override;

    // output
    TTree* fTree = 0;

    enum SynchPoint {
      kGoodVertex,
      kCleanJets,
      kLeadingJet110,
      kDeltaPhiJ1J2,
      kMET200,
      kNAK4Jets,
      kEMuVeto,
      kTauVeto,
      kPhotonVeto,   
      nSynchPoints
    };
    Bool_t fSynchPass[nSynchPoints] = {};

    // input
    TString fVerticesName = "Vertices";
    TString fMetName = "PFMet";
    TString fJetsName = "AKt4PFJets";
    TString fElectronsName = "Electrons";
    TString fMuonsName = "Muons";
    TString fTausName = "HPSTaus";
    TString fPhotonsName = "Photons";

    ClassDef(PreRun2SynchExercise, 0)
  };

}

#endif

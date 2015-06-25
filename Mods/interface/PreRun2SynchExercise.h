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

    char const* SetVerticesName() const { return fVerticesName; }
    char const* SetMetName() const { return fMetName; }
    char const* SetJetsName() const { return fJetsName; }
    char const* SetElectronsName() const { return fElectronsName; }
    char const* SetMuonsName() const { return fMuonsName; }
    char const* SetTausName() const { return fTausName; }
    char const* SetPhotonsName() const { return fPhotonsName; }

    void SetVerticesName(char const* n) { fVerticesName = n; }
    void SetMetName(char const* n) { fMetName = n; }
    void SetJetsName(char const* n) { fJetsName = n; }
    void SetElectronsName(char const* n) { fElectronsName = n; }
    void SetMuonsName(char const* n) { fMuonsName = n; }
    void SetTausName(char const* n) { fTausName = n; }
    void SetPhotonsName(char const* n) { fPhotonsName = n; }

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

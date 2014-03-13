//--------------------------------------------------------------------------------------------------
// BoostedVTreeWriter
//
// This module writes a tree about jets corresponding to Vector bosons. See documentation here:
//
//
// Authors: C.Paus
//--------------------------------------------------------------------------------------------------
#ifndef MITMONOJET_MODS_BOOSTEDVTREEWRITER_H
#define MITMONOJET_MODS_BOOSTEDVTREEWRITER_H

#include <TFile.h>
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
#include "MitAna/DataTree/interface/PFTauCol.h"
#include "MitAna/DataTree/interface/PFCandidateCol.h"
#include "MitMonoJet/Core/MitGPBoostedVTree.h"

namespace mithep
{

  class BoostedVTreeWriter : public BaseMod
  {
  public:
    typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > LorentzVector;

    BoostedVTreeWriter(const char *name  = "BoostedVTreeWriter",
		       const char *title = "Vector Boson Tagging module");

    ~BoostedVTreeWriter();

    void                          SetIsData(Bool_t b)               { fIsData = b; }
    void                          SetMcPartsName(const char *n)     { fMcPartsName = n; }
    void                          SetTriggerObjsName(const char *n) { fTriggerObjsName = n; }
    void                          SetJetsName(const char *n)        { fJetsName = n; }
    void                          SetJetsFromBranch(bool b)         { fJetsFromBranch = b; }
    void                          SetPFCandidatesName(const char *n){ fPFCandidatesName = n; }
    void                          SetPFCandidatesFromBranch(bool b) { fPFCandidatesFromBranch = b; }
    void                          SetPhotonsName(const char *n)     { fPhotonsName= n; }
    void                          SetPhotonsFromBranch(bool b)      { fPhotonsFromBranch = b; }
    void                          SetPFTausName(const char *n)      { fPFTausName = n; }
    void                          SetPFTausFromBranch(bool b)       { fPFTausFromBranch = b; }
    void                          SetLeptonsName(const char *n)     { fLeptonsName  = n; } 
    void                          SetPFNoPileUpName(const char *n)  { fPFNoPileUpName  = n; } 
    void                          SetPFPileUpName(const char *n)    { fPFPileUpName  = n; }
    void                          SetOutputName(const char *n)      { fOutputName = n; }
    void                          SetPruning(Int_t n)               { fPrune = n; }

    void                          SetHistNPtBins(Int_t n)           { fHistNPtBins = n; }
    void                          SetHistNEtaBins(Int_t n)          { fHistNEtaBins = n; }
    void                          SetHistMinPt(Double_t b)          { fHistMinPt = b; }
    void                          SetHistMaxPt(Double_t b)          { fHistMaxPt = b; }
    void                          SetHistMinEta(Double_t b)         { fHistMinEta = b; }
    void                          SetHistMaxEta(Double_t b)         { fHistMaxEta = b; }
    void                          SetHistTau1Bins(Int_t n)          { fHistTau1Bins = n; }
    void                          SetHistTau2Bins(Int_t n)          { fHistTau2Bins = n; }
    void                          SetHistTau3Bins(Int_t n)          { fHistTau3Bins = n; }
    void                          SetHistT2ovrT1Bins(Double_t b)    { fHistT2ovrT1Bins = b; }
    void                          SetHistT3ovrT2Bins(Double_t b)    { fHistT3ovrT2Bins = b; }
    void                          SetHistMinTau1(Double_t b)        { fHistMinTau1 = b; }
    void                          SetHistMinTau2(Double_t b)        { fHistMinTau2 = b; }
    void                          SetHistMinTau3(Double_t b)        { fHistMinTau3 = b; }
    void                          SetHistMinT2ovrT1(Double_t b)     { fHistMinT2ovrT1 = b; }
    void                          SetHistMinT3ovrT2(Double_t b)     { fHistMinT3ovrT2 = b; }
    void                          SetHistMaxTau1(Double_t b)        { fHistMaxTau1 = b; }
    void                          SetHistMaxTau2(Double_t b)        { fHistMaxTau2 = b; }
    void                          SetHistMaxTau3(Double_t b)        { fHistMaxTau3 = b; }
    void                          SetHistMaxT2ovrT1(Double_t b)     { fHistMaxT2ovrT1 = b; }
    void                          SetHistMaxT3ovrT2(Double_t b)     { fHistMaxT3ovrT2 = b; }

  protected:
    void                          Process();
    void                          SlaveBegin();
    void                          SlaveTerminate();
  private:
    void                          ProcessMc();             // Deals with the MC separately
    void                          GetJetTriggerObjs();     // Fills the jet trigger objs
    Double_t                      MinTriggerDeltaR(LorentzVector jet);
    float                         GetTau(fastjet::PseudoJet &iJet,int iN, float iKappa);
    bool                          IsTightMuon(const Muon *muon);

    Bool_t                        fIsData;                 //is this data or MC?
    TString                       fMcPartsName;            //(i) name of MC particles
    const MCParticleCol          *fMcParts;	           //MC particle coll
    TString                       fTriggerObjsName;        //(i) name of trigger objects
    const TriggerObjectCol       *fTrigObjs;               //trigger objects coll handle
    TString                       fJetsName;               //(i) name of jets used to make trigger
    Bool_t                        fJetsFromBranch;         //are jets from Branch?
    const JetCol                 *fJets;                   //jets used to make the trigger
    TString                       fPFCandidatesName;       //(i) name of PF candidates coll
    Bool_t                        fPFCandidatesFromBranch; //are PF candidates from Branch?
    const PFCandidateCol         *fPFCandidates;           //particle flow candidates coll handle
    TString                       fPhotonsName;            //(i) name of photon coll
    Bool_t                        fPhotonsFromBranch;      //are photons from Branch?
    const PhotonCol              *fPhotons;                //photon coll handle
    TString                       fPFTausName;             //(i) name of PF tau coll
    Bool_t                        fPFTausFromBranch;       //are PF taus from Branch?
    const PFTauCol               *fPFTaus;                 //PF tau coll handle

    TString                       fLeptonsName;            //(i) name of the merged leptons
    TString                       fPFNoPileUpName;         //(i) name for PF no pileup candidates
    TString                       fPFPileUpName;           //(i) name for PF pileup candidates

    std::vector<const TriggerObject *>
                                  fJetTriggerObjs;         //jet trigger objects

    // Objects from fastjet we want to use
    double                        fConeSize;
    int                           fPrune;                  //apply pruning: 0-no, 1-standard CMS
    fastjet::Pruner              *fPruner;
    fastjet::JetDefinition       *fCAJetDef;
    fastjet::GhostedAreaSpec     *fActiveArea;
    fastjet::AreaDefinition      *fAreaDefinition;

    // Output histograms
    TH1D                         *fPFCandidatesPt;         //particle flow pt  distribution
    TH1D                         *fPFCandidatesEta;        //particle flow eta distribution
    TH1D                         *fCAJetPt;                //fastjet output pt distribution
    TH1D                         *fCAJetEta;               //fastjet output eta distribution
    TH1D                         *fCATau1;                 //plot of Tau 1
    TH1D                         *fCATau2;                 //plot of Tau 2
    TH1D                         *fCATau3;                 //plot of Tau 3
    TH1D                         *fCAT2ovrT1;              //plot of Tau 2 over Tau 1
    TH1D                         *fCAT3ovrT2;              //plot of Tau 3 over Tau 2

    Int_t                         fNAnalyzed;              //number of events analyzed

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
    TString                       fOutputName;             //(o) name of ntuple output
    TFile	                 *fOutputFile;
    MitGPBoostedVTree             fMitGPTree;

    ClassDef(BoostedVTreeWriter, 0) // Boosted Vector boson tree writer
  };
}
#endif

//--------------------------------------------------------------------------------------------------
// MitGPBoostedVTree
//
// This is the definition that we use to write a tree with information relevant to boosted Vector
// bosons.
//
// Authors: C.Freeman, O.Onuoha, L.DiMatteo, C.Paus
//--------------------------------------------------------------------------------------------------
#ifndef MitGPBoostedVTree_H
#define MitGPBoostedVTree_H

#include <set>
#include <vector>
#include <string>
#include <utility>

#include <TFile.h>
#include <TTree.h>
#include <TNtuple.h>
#include <TError.h>

#include "Math/VectorUtil.h"
#include "Math/LorentzVector.h"

class MitGPBoostedVTree
{
public:
  // float doesn't have dictionary by default, so use double
  typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > LorentzVector;

  /// variables
  unsigned int           nPFCandidates_;
  float                  jet1Eta_;
  float                  jet1Phi_;
  float                  jet1Pt_;
  float                  jet1Tau1_;
  float                  jet1Tau2_;
  float                  jet1Tau3_;
  float                  jet1R_;
  float                  jet1M_;
  float                  jet2Eta_;
  float                  jet2Phi_;
  float                  jet2Pt_;
  float                  jet2Tau1_;
  float                  jet2Tau2_;
  float                  jet2Tau3_;
  float                  jet2R_;
  float                  jet2M_;
  float                  eta_;
  float                  phi_;
  float                  pt_;
  int                    numJets_;
  LorentzVector          jet1_;
  LorentzVector          jet2_;

  // this is the main element
  TTree                 *tree_;
  TFile                 *file_;

  // hold the names of variables to facilitate things (filled during Init)
  std::vector<std::string> variables_;

  // default constructor
  MitGPBoostedVTree(): jetPtr1_(&jet1_), jetPtr2_(&jet2_){}
  // default destructor
  ~MitGPBoostedVTree()  { if (file_) file_->Close(); };

  /// initialize varibles and fill list of available variables
  void InitVariables();

  /// load a MitGPBoostedVTree
  void LoadTree(const char* file, int type = -1){
    // to load three different ntuples in the same job DMTree0/1
    // type == 0/1 if all variables was added
    // type = -1 (default) if a minimum set of variables was added with tree as name
    file_ = TFile::Open(file);
    assert(file_);
    if     (type == 0)
      tree_ = dynamic_cast<TTree*>(file_->FindObjectAny("MJetTree"));
    else if(type == 1)
      tree_ = dynamic_cast<TTree*>(file_->FindObjectAny("MJetTreeDiLepton"));
    else if(type == 2)
      tree_ = dynamic_cast<TTree*>(file_->FindObjectAny("MJetTreeWlnu"));
    else
      tree_ = dynamic_cast<TTree*>(file_->FindObjectAny("tree"));
    assert(tree_);
  }

  // create a MitGPBoostedVTree
  void CreateTree(int type = -1)
  {
    assert(type==type);  // just to suppress warnings

    // to create three different ntuples in the same job DMTree0/1
    // type == 0/1 add all variables
    // type = -1 (default) add a minimum set of variables with tree as name

    if     (type == 0)
      tree_ = new TTree("MJetTree","Smurf ntuples");
    else if(type == 1)
      tree_ = new TTree("MJetTreeDiLepton","Smurf ntuples");
    else if(type == 2)
      tree_ = new TTree("MJetTreeWlnu","Smurf ntuples");
    else
      tree_ = new TTree("tree","Smurf ntuples");

    file_ = 0;

    InitVariables();

    //book the branches
    tree_->Branch("nPFCandidates", &nPFCandidates_,"nPFCandidates/i");
    tree_->Branch("jet1Eta"      , &jet1Eta_,      "jet1Eta/F");
    tree_->Branch("jet1Phi"      , &jet1Phi_,      "jet1Phi/F");
    tree_->Branch("jet1Pt"       , &jet1Pt_,       "jet1Pt/F");
    tree_->Branch("Eta"          , &eta_,          "Eta/F");
    tree_->Branch("Phi"          , &phi_,          "Phi/F");
    tree_->Branch("Pt"           , &pt_,           "Pt/F");
    tree_->Branch("numJets"      , &numJets_,      "numJets/i");
    tree_->Branch("jet1R"        , &jet1R_,        "jet1R/F");
    tree_->Branch("jet1M"        , &jet1M_,        "jet1M/F");
    tree_->Branch("jet1Tau1"     , &jet1Tau1_,     "jet1Tau1/F");
    tree_->Branch("jet1Tau2"     , &jet1Tau2_,     "jet1Tau2/F");
    tree_->Branch("jet1Tau3"     , &jet1Tau3_,     "jet1Tau3/F");
    tree_->Branch("jet2Eta"      , &jet2Eta_,      "jet2Eta/F");
    tree_->Branch("jet2Phi"      , &jet2Phi_,      "jet2Phi/F");
    tree_->Branch("jet2Pt"       , &jet2Pt_,       "jet2Pt/F");
    tree_->Branch("jet2R"        , &jet2R_,        "jet2R/F");
    tree_->Branch("jet2M"        , &jet2M_,        "jet2M/F");
    tree_->Branch("jet2Tau1"     , &jet2Tau1_,     "jet2Tau1/F");
    tree_->Branch("jet2Tau2"     , &jet2Tau2_,     "jet2Tau2/F");
    tree_->Branch("jet2Tau3"     , &jet2Tau3_,     "jet2Tau3/F");
    tree_->Branch("jet1"         , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >", &jetPtr1_);
    tree_->Branch("jet2"         , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >", &jetPtr2_);
  }

  // initialze a MitGPBoostedVTree
  void InitTree(int type = -1)
  {
    assert(type==type); // just to suppress warnings
    assert(tree_);

    // don't forget to set pointers to zero before you set address
    // or you will fully appreciate that "ROOT sucks" :)
    InitVariables();

    //Set branch address
    Int_t currentState = gErrorIgnoreLevel;
    // gErrorIgnoreLevel = kError;
    gErrorIgnoreLevel = kBreak;
    tree_->SetBranchAddress("nPFCandidates", &nPFCandidates_);
    tree_->SetBranchAddress("jet1Eta"      , &jet1Eta_);
    tree_->SetBranchAddress("jet1Phi"      , &jet1Phi_);
    tree_->SetBranchAddress("jet1Pt"       , &jet1Pt_);
    tree_->SetBranchAddress("Eta"          , &eta_);
    tree_->SetBranchAddress("Phi"          , &phi_);
    tree_->SetBranchAddress("Pt"           , &pt_);
    tree_->SetBranchAddress("numJets"      , &numJets_);
    tree_->SetBranchAddress("jet1R"        , &jet1R_);
    tree_->SetBranchAddress("jet1M"        , &jet1M_);
    tree_->SetBranchAddress("jet1Tau1"     , &jet1Tau1_);
    tree_->SetBranchAddress("jet1Tau2"     , &jet1Tau2_);
    tree_->SetBranchAddress("jet1Tau3"     , &jet1Tau3_);
    tree_->SetBranchAddress("jet2Eta"      , &jet2Eta_);
    tree_->SetBranchAddress("jet2Phi"      , &jet2Phi_);
    tree_->SetBranchAddress("jet2Pt"       , &jet2Pt_);
    tree_->SetBranchAddress("jet2R"        , &jet2R_);
    tree_->SetBranchAddress("jet2M"        , &jet2M_);
    tree_->SetBranchAddress("jet2Tau1"     , &jet2Tau1_);
    tree_->SetBranchAddress("jet2Tau2"     , &jet2Tau2_);
    tree_->SetBranchAddress("jet2Tau3"     , &jet2Tau3_);
    tree_->SetBranchAddress("jet1"         , &jetPtr1_);
    tree_->SetBranchAddress("jet2"         , &jetPtr2_);

    gErrorIgnoreLevel = currentState;
  }

private:
  LorentzVector* jetPtr1_;
  LorentzVector* jetPtr2_;
};

inline void
MitGPBoostedVTree::InitVariables()
{
  // inizialize variables
  nPFCandidates_  = 0;
  numJets_        = 0;
  jet1Eta_        = 0;
  jet1Phi_        = 0;
  jet1Pt_         = 0;
  jet1Tau1_       = 0;
  jet1Tau2_       = 0;
  jet1Tau3_       = 0;
  jet1R_          = 0;
  jet1M_          = 0;
  jet2Eta_        = 0;
  jet2Phi_        = 0;
  jet2Pt_         = 0;
  jet2Tau1_       = 0;
  jet2Tau2_       = 0;
  jet2Tau3_       = 0;
  jet2R_          = 0;
  jet2M_          = 0;
  jet1_           = LorentzVector();
  jet2_           = LorentzVector();
}

#endif

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

  // variables
  unsigned int           nGenParts_;
  int                    numGenJets_;

  LorentzVector          genJet1_;
  float                  genJet1NParts_;
  float                  genJet1Eta_;
  float                  genJet1Phi_;
  float                  genJet1Pt_;
  float                  genJet1Tau1_;
  float                  genJet1Tau2_;
  float                  genJet1Tau3_;
  float                  genJet1R_;
  float                  genJet1M_;
  float                  genJet1MinTrigDr_;

  LorentzVector          genJet2_;
  float                  genJet2NParts_;
  float                  genJet2Eta_;
  float                  genJet2Phi_;
  float                  genJet2Pt_;
  float                  genJet2Tau1_;
  float                  genJet2Tau2_;
  float                  genJet2Tau3_;
  float                  genJet2R_;
  float                  genJet2M_;
  float                  genJet2MinTrigDr_;

  unsigned int           event_;
  unsigned int           run_;
  unsigned int           lumi_;
  unsigned int           trigger_;
  unsigned int           nParts_;
  int                    numJets_;

  LorentzVector          jet1_;
  float                  jet1NParts_;
  float                  jet1Eta_;
  float                  jet1Phi_;
  float                  jet1Pt_;
  float                  jet1Tau1_;
  float                  jet1Tau2_;
  float                  jet1Tau3_;
  float                  jet1R_;
  float                  jet1M_;
  float                  jet1MinTrigDr_;
  float                  jet1QGTag_;

  LorentzVector          jet2_;
  float                  jet2NParts_;
  float                  jet2Eta_;
  float                  jet2Phi_;
  float                  jet2Pt_;
  float                  jet2Tau1_;
  float                  jet2Tau2_;
  float                  jet2Tau3_;
  float                  jet2R_;
  float                  jet2M_;
  float                  jet2MinTrigDr_;
  float                  jet2QGTag_;

  unsigned int           nlep_;
  LorentzVector          lep1_;
  int                    lid1_;
  int		         lep1IsTightMuon_;
  int		         lep1IsIsolated_;
  float                  lep1PtErr_;
  LorentzVector          lep2_;
  int                    lid2_;
  int		         lep2IsTightMuon_;
  int		         lep2IsIsolated_;
  float                  lep2PtErr_;
  LorentzVector          lep3_;
  int                    lid3_;
  int		         lep3IsTightMuon_;
		         
  unsigned int           ntaus_;
  LorentzVector          tau1_;
  LorentzVector          tau2_;
		         
  unsigned int           nphotons_;
  LorentzVector          pho1_;

  // this is the main element
  TTree                 *tree_;
  TFile                 *file_;

  // hold the names of variables to facilitate things (filled during Init)
  std::vector<std::string> variables_;

  // default constructor
  MitGPBoostedVTree() :
    jetPtr1_(&jet1_), jetPtr2_(&jet2_),
    genJetPtr1_(&genJet1_), genJetPtr2_(&genJet2_),
    lepPtr1_(&lep1_), lepPtr2_(&lep2_), lepPtr3_(&lep3_),
    tauPtr1_(&tau1_), tauPtr2_(&tau2_),
    phoPtr1_(&pho1_)
  {}

  // default destructor
  ~MitGPBoostedVTree() { if (file_) file_->Close(); };

  /// initialize varibles and fill list of available variables
  void InitVariables();

  // load a MitGPBoostedVTree
  void LoadTree(const char* file, int type = -1) {
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
    tree_->Branch("event"           , &event_,           "event/i");
    tree_->Branch("run"             , &run_,             "run/i");
    tree_->Branch("lumi"            , &lumi_,            "lumi/i");
    tree_->Branch("trigger"         , &trigger_,         "trigger/i");

    tree_->Branch("nGenParts"       , &nGenParts_,       "nGenParts/i");
    tree_->Branch("numGenJets"      , &numGenJets_,      "numGenJets/i");

    tree_->Branch("genJet1"         , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >", &genJetPtr1_);
    tree_->Branch("genJet1NParts"   , &genJet1NParts_,   "genJet1NParts/I");
    tree_->Branch("genJet1Eta"      , &genJet1Eta_,      "genJet1Eta/F");
    tree_->Branch("genJet1Phi"      , &genJet1Phi_,      "genJet1Phi/F");
    tree_->Branch("genJet1Pt"       , &genJet1Pt_,       "genJet1Pt/F");
    tree_->Branch("genJet1R"        , &genJet1R_,        "genJet1R/F");
    tree_->Branch("genJet1M"        , &genJet1M_,        "genJet1M/F");
    tree_->Branch("genJet1Tau1"     , &genJet1Tau1_,     "genJet1Tau1/F");
    tree_->Branch("genJet1Tau2"     , &genJet1Tau2_,     "genJet1Tau2/F");
    tree_->Branch("genJet1Tau3"     , &genJet1Tau3_,     "genJet1Tau3/F");
    tree_->Branch("genJet1MinTrigDr", &genJet1MinTrigDr_,"genJet1MinTrigDr/F");

    tree_->Branch("genJet2"         , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >", &genJetPtr2_);
    tree_->Branch("genJet2NParts"   , &genJet2NParts_,   "genJet2NParts/I");
    tree_->Branch("genJet2Eta"      , &genJet2Eta_,      "genJet2Eta/F");
    tree_->Branch("genJet2Phi"      , &genJet2Phi_,      "genJet2Phi/F");
    tree_->Branch("genJet2Pt"       , &genJet2Pt_,       "genJet2Pt/F");
    tree_->Branch("genJet2R"        , &genJet2R_,        "genJet2R/F");
    tree_->Branch("genJet2M"        , &genJet2M_,        "genJet2M/F");
    tree_->Branch("genJet2Tau1"     , &genJet2Tau1_,     "genJet2Tau1/F");
    tree_->Branch("genJet2Tau2"     , &genJet2Tau2_,     "genJet2Tau2/F");
    tree_->Branch("genJet2Tau3"     , &genJet2Tau3_,     "genJet2Tau3/F");
    tree_->Branch("genJet2MinTrigDr", &genJet2MinTrigDr_,"genJet2MinTrigDr/F");

    tree_->Branch("nParts"          , &nParts_,          "nParts/i");
    tree_->Branch("numJets"         , &numJets_,         "numJets/i");

    tree_->Branch("jet1"            , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >", &jetPtr1_);
    tree_->Branch("jet1NParts"      , &jet1NParts_,      "jet1NParts/I");
    tree_->Branch("jet1Eta"         , &jet1Eta_,         "jet1Eta/F");
    tree_->Branch("jet1Phi"         , &jet1Phi_,         "jet1Phi/F");
    tree_->Branch("jet1Pt"          , &jet1Pt_,          "jet1Pt/F");
    tree_->Branch("jet1R"           , &jet1R_,           "jet1R/F");
    tree_->Branch("jet1M"           , &jet1M_,           "jet1M/F");
    tree_->Branch("jet1Tau1"        , &jet1Tau1_,        "jet1Tau1/F");
    tree_->Branch("jet1Tau2"        , &jet1Tau2_,        "jet1Tau2/F");
    tree_->Branch("jet1Tau3"        , &jet1Tau3_,        "jet1Tau3/F");
    tree_->Branch("jet1MinTrigDr"   , &jet1MinTrigDr_,   "jet1MinTrigDr/F");
    tree_->Branch("jet1QGTag"       , &jet1QGTag_,       "jet1QGTag/F");

    tree_->Branch("jet2"            , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >", &jetPtr2_);
    tree_->Branch("jet2NParts"      , &jet2NParts_,      "jet2NParts/I");
    tree_->Branch("jet2Eta"         , &jet2Eta_,         "jet2Eta/F");
    tree_->Branch("jet2Phi"         , &jet2Phi_,         "jet2Phi/F");
    tree_->Branch("jet2Pt"          , &jet2Pt_,          "jet2Pt/F");
    tree_->Branch("jet2R"           , &jet2R_,           "jet2R/F");
    tree_->Branch("jet2M"           , &jet2M_,           "jet2M/F");
    tree_->Branch("jet2Tau1"        , &jet2Tau1_,        "jet2Tau1/F");
    tree_->Branch("jet2Tau2"        , &jet2Tau2_,        "jet2Tau2/F");
    tree_->Branch("jet2Tau3"        , &jet2Tau3_,        "jet2Tau3/F");
    tree_->Branch("jet2MinTrigDr"   , &jet2MinTrigDr_,   "jet2MinTrigDr/F");
    tree_->Branch("jet2QGTag"       , &jet2QGTag_,       "jet2QGTag/F");

    tree_->Branch("nlep"            , &nlep_,            "nlep/i");
    tree_->Branch("lep1"            , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >", &lepPtr1_);
    tree_->Branch("lid1"            , &lid1_,            "lid1/I");
    tree_->Branch("lep1IsTightMuon" , &lep1IsTightMuon_, "lep1IsTightMuon/I");
    tree_->Branch("lep1IsIsolated"  , &lep1IsIsolated_,  "lep1IsIsolated/I");
    tree_->Branch("lep1PtErr"       , &lep1PtErr_,       "lep1PtErr/F");
    tree_->Branch("lep2"            , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >", &lepPtr2_);
    tree_->Branch("lid2"            , &lid2_,            "lid2/I");
    tree_->Branch("lep2IsTightMuon" , &lep2IsTightMuon_, "lep2IsTightMuon/I");
    tree_->Branch("lep2IsIsolated"  , &lep2IsIsolated_,  "lep2IsIsolated/I");
    tree_->Branch("lep2PtErr"       , &lep2PtErr_,       "lep2PtErr/F");
    tree_->Branch("lep3"            , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >", &lepPtr3_);
    tree_->Branch("lid3"            , &lid3_,            "lid3/I");
    tree_->Branch("lep3IsTightMuon" , &lep3IsTightMuon_, "lep3IsTightMuon/I");

    tree_->Branch("ntaus"           , &ntaus_,           "ntaus/i");
    tree_->Branch("tau1"            , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >", &tauPtr1_);
    tree_->Branch("tau2"            , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >", &tauPtr2_);

    tree_->Branch("nphotons"        , &nphotons_,        "nphotons/i");
    tree_->Branch("pho1"            , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >", &phoPtr1_);
  }

  // initialze a MitGPBoostedVTree
  void InitTree(int type = -1)
  {
    assert(type==type); // just to suppress warnings
    assert(tree_);      // make sure the tree exists

    // Don't forget to set pointers to zero before you set address or you will fully appreciate that
    // "ROOT sucks" :)
    InitVariables();

    //Set branch address
    Int_t currentState = gErrorIgnoreLevel;
    gErrorIgnoreLevel = kBreak;

    tree_->SetBranchAddress("event"           , &event_);
    tree_->SetBranchAddress("run"             , &run_);
    tree_->SetBranchAddress("lumi"            , &lumi_);
    tree_->SetBranchAddress("trigger"         , &trigger_);

    tree_->SetBranchAddress("nGenParts"       , &nGenParts_);
    tree_->SetBranchAddress("numGenJets"      , &numGenJets_);

    tree_->SetBranchAddress("genJet1"         , &genJetPtr1_);
    tree_->SetBranchAddress("genJet1NParts"   , &genJet1NParts_);
    tree_->SetBranchAddress("genJet1Eta"      , &genJet1Eta_);
    tree_->SetBranchAddress("genJet1Phi"      , &genJet1Phi_);
    tree_->SetBranchAddress("genJet1Pt"       , &genJet1Pt_);
    tree_->SetBranchAddress("genJet1R"        , &genJet1R_);
    tree_->SetBranchAddress("genJet1M"        , &genJet1M_);
    tree_->SetBranchAddress("genJet1Tau1"     , &genJet1Tau1_);
    tree_->SetBranchAddress("genJet1Tau2"     , &genJet1Tau2_);
    tree_->SetBranchAddress("genJet1Tau3"     , &genJet1Tau3_);
    tree_->SetBranchAddress("genJet1MinTrigDr", &genJet1MinTrigDr_);

    tree_->SetBranchAddress("genJet2"         , &genJetPtr2_);
    tree_->SetBranchAddress("genJet2NParts"   , &genJet2NParts_);
    tree_->SetBranchAddress("genJet2Eta"      , &genJet2Eta_);
    tree_->SetBranchAddress("genJet2Phi"      , &genJet2Phi_);
    tree_->SetBranchAddress("genJet2Pt"       , &genJet2Pt_);
    tree_->SetBranchAddress("genJet2R"        , &genJet2R_);
    tree_->SetBranchAddress("genJet2M"        , &genJet2M_);
    tree_->SetBranchAddress("genJet2Tau1"     , &genJet2Tau1_);
    tree_->SetBranchAddress("genJet2Tau2"     , &genJet2Tau2_);
    tree_->SetBranchAddress("genJet2Tau3"     , &genJet2Tau3_);
    tree_->SetBranchAddress("genJet2MinTrigDr", &genJet2MinTrigDr_);

    tree_->SetBranchAddress("nParts"          , &nParts_);
    tree_->SetBranchAddress("numJets"         , &numJets_);

    tree_->SetBranchAddress("jet1"            , &jetPtr1_);
    tree_->SetBranchAddress("jet1NParts"      , &jet1NParts_);
    tree_->SetBranchAddress("jet1Eta"         , &jet1Eta_);
    tree_->SetBranchAddress("jet1Phi"         , &jet1Phi_);
    tree_->SetBranchAddress("jet1Pt"          , &jet1Pt_);
    tree_->SetBranchAddress("jet1R"           , &jet1R_);
    tree_->SetBranchAddress("jet1M"           , &jet1M_);
    tree_->SetBranchAddress("jet1Tau1"        , &jet1Tau1_);
    tree_->SetBranchAddress("jet1Tau2"        , &jet1Tau2_);
    tree_->SetBranchAddress("jet1Tau3"        , &jet1Tau3_);
    tree_->SetBranchAddress("jet1MinTrigDr"   , &jet1MinTrigDr_);
    tree_->SetBranchAddress("jet1QGTag"       , &jet1QGTag_);

    tree_->SetBranchAddress("jet2"            , &jetPtr2_);
    tree_->SetBranchAddress("jet2NParts"      , &jet2NParts_);
    tree_->SetBranchAddress("jet2Eta"         , &jet2Eta_);
    tree_->SetBranchAddress("jet2Phi"         , &jet2Phi_);
    tree_->SetBranchAddress("jet2Pt"          , &jet2Pt_);
    tree_->SetBranchAddress("jet2R"           , &jet2R_);
    tree_->SetBranchAddress("jet2M"           , &jet2M_);
    tree_->SetBranchAddress("jet2Tau1"        , &jet2Tau1_);
    tree_->SetBranchAddress("jet2Tau2"        , &jet2Tau2_);
    tree_->SetBranchAddress("jet2Tau3"        , &jet2Tau3_);
    tree_->SetBranchAddress("jet2MinTrigDr"   , &jet2MinTrigDr_);
    tree_->SetBranchAddress("jet2QGTag"       , &jet2QGTag_);

    tree_->SetBranchAddress("nlep"            , &nlep_);
    tree_->SetBranchAddress("lep1"            , &lepPtr1_);
    tree_->SetBranchAddress("lid1"            , &lid1_);
    tree_->SetBranchAddress("lep1IsTightMuon" , &lep1IsTightMuon_);
    tree_->SetBranchAddress("lep1IsIsolated"  , &lep1IsIsolated_);
    tree_->SetBranchAddress("lep1PtErr"       , &lep1PtErr_);
    tree_->SetBranchAddress("lep2"            , &lepPtr2_);
    tree_->SetBranchAddress("lid2"            , &lid2_);
    tree_->SetBranchAddress("lep2IsTightMuon" , &lep2IsTightMuon_);
    tree_->SetBranchAddress("lep2IsIsolated"  , &lep2IsIsolated_);
    tree_->SetBranchAddress("lep2PtErr"       , &lep2PtErr_);
    tree_->SetBranchAddress("lep3"            , &lepPtr3_);
    tree_->SetBranchAddress("lid3"            , &lid3_);
    tree_->SetBranchAddress("lep3IsTightMuon" , &lep3IsTightMuon_);

    tree_->SetBranchAddress("ntaus"           , &ntaus_);
    tree_->SetBranchAddress("tau1"            , &tauPtr1_);
    tree_->SetBranchAddress("tau2"            , &tauPtr2_);

    tree_->SetBranchAddress("nphotons"        , &nphotons_);
    tree_->SetBranchAddress("pho1"            , &phoPtr1_);

    gErrorIgnoreLevel = currentState;
  }

private:
  LorentzVector* jetPtr1_;
  LorentzVector* jetPtr2_;
  LorentzVector* genJetPtr1_;
  LorentzVector* genJetPtr2_;
  LorentzVector* lepPtr1_;
  LorentzVector* lepPtr2_;
  LorentzVector* lepPtr3_;
  LorentzVector* tauPtr1_;
  LorentzVector* tauPtr2_;
  LorentzVector* phoPtr1_;
};

inline void
MitGPBoostedVTree::InitVariables()
{
  // inizialize variables
  event_            = 0;
  run_              = 0;
  lumi_             = 0;
  trigger_          = 0;

  nGenParts_        = 0;
  numGenJets_       = 0;

  genJet1_          = LorentzVector();
  genJet1Eta_       = 0;
  genJet1NParts_    = 0;
  genJet1Phi_       = 0;
  genJet1Pt_        = 0;
  genJet1Tau1_      = 0;
  genJet1Tau2_      = 0;
  genJet1Tau3_      = 0;
  genJet1R_         = 0;
  genJet1M_         = 0;
  genJet1MinTrigDr_ = 999;

  genJet2_          = LorentzVector();
  genJet2NParts_    = 0;
  genJet2Eta_       = 0;
  genJet2Phi_       = 0;
  genJet2Pt_        = 0;
  genJet2Tau1_      = 0;
  genJet2Tau2_      = 0;
  genJet2Tau3_      = 0;
  genJet2R_         = 0;
  genJet2M_         = 0;
  genJet2MinTrigDr_ = 999;

  nParts_           = 0;
  numJets_          = 0;

  jet1_             = LorentzVector();
  jet1NParts_       = 0;
  jet1Eta_          = 0;
  jet1Phi_          = 0;
  jet1Pt_           = 0;
  jet1Tau1_         = 0;
  jet1Tau2_         = 0;
  jet1Tau3_         = 0;
  jet1R_            = 0;
  jet1M_            = 0;
  jet1MinTrigDr_    = 999;
  jet1QGTag_        = 0;

  jet2_             = LorentzVector();
  jet2NParts_       = 0;
  jet2Eta_          = 0;
  jet2Phi_          = 0;
  jet2Pt_           = 0;
  jet2Tau1_         = 0;
  jet2Tau2_         = 0;
  jet2Tau3_         = 0;
  jet2R_            = 0;
  jet2M_            = 0;
  jet2MinTrigDr_    = 999;
  jet2QGTag_        = 0;

  nlep_             = 0;
  lep1_       	    = LorentzVector();
  lid1_             = 0;
  lep1IsTightMuon_  = 0;
  lep1IsIsolated_   = 0;
  lep1PtErr_        = 999;
  lep2_       	    = LorentzVector();
  lid2_             = 0;
  lep2IsTightMuon_  = 0;
  lep2IsIsolated_   = 0;
  lep2PtErr_        = 999;
  lep3_       	    = LorentzVector();
  lid3_             = 0;
  lep3IsTightMuon_  = 0;

  ntaus_            = 0;
  tau1_       	    = LorentzVector();
  tau2_       	    = LorentzVector();

  nphotons_         = 0;
  pho1_       	    = LorentzVector();
}

#endif

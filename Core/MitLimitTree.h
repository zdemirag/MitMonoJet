#ifndef MitLimitTree_H
#define MitLimitTree_H

#include <set>
#include <vector>
#include <string>
#include <utility>

#include "TFile.h"
#include "TTree.h"
#include "TError.h"

class MitLimitTree {
 public:

  /// variables
  Float_t          mvamet_;
  Float_t          mvametphi_;
  Float_t          jet1pt_;
  Float_t          genjetpt_;
  Float_t          genVpt_;
  Float_t          dmpt_;
  Float_t          weight_;
  
 public:
  /// this is the main element
  TFile *f_;
  TTree *tree_;

  /// default constructor  
  MitLimitTree() {}
  /// default destructor
  ~MitLimitTree(){
    if (f_) f_->Close();  
  };

  /// initialize varibles and fill list of available variables
  void InitVariables();

  /// load a MitLimitTree
  void LoadTree(const char* file, const char* treeName){
    f_ = TFile::Open(file);
    assert(f_);
    tree_ = dynamic_cast<TTree*>(f_->FindObjectAny(treeName));
    assert(tree_);
  }

  /// create a MitLimitTree
  void CreateTree(const char* treeName){
    tree_ = new TTree(treeName,"Limit tree");
    f_ = 0;
    InitVariables();
    
    //book the branches    
    tree_->Branch("mvamet"         , &mvamet_          ,   "mvamet/F");
    tree_->Branch("mvametphi"      , &mvametphi_       ,   "mvametphi/F");
    tree_->Branch("jet1pt"         , &jet1pt_          ,   "jet1pt/F");
    tree_->Branch("genjetpt"       , &genjetpt_        ,   "genjetpt/F");
    tree_->Branch("genVpt"         , &genVpt_          ,   "genVpt/F");
    tree_->Branch("dmpt"           , &dmpt_            ,   "dmpt/F");
    tree_->Branch("weight"         , &weight_          ,   "weight/F");
  }

  /// create a MitLimitTree
  void DestroyTree(){
    tree_->Delete();
  }

  // initialze a MitLimitTree
  void InitTree(){
    assert(tree_);
    // don't forget to set pointers to zero before you set address
    // or you will fully appreciate that "ROOT sucks" :)
    InitVariables();
    //Set branch address
    Int_t currentState = gErrorIgnoreLevel;
    // gErrorIgnoreLevel = kError;
    gErrorIgnoreLevel = kBreak;
    tree_->SetBranchAddress("mvamet"         , &mvamet_);
    tree_->SetBranchAddress("mvametphi"      , &mvametphi_);
    tree_->SetBranchAddress("jet1pt"         , &jet1pt_);
    tree_->SetBranchAddress("genjetpt"       , &genjetpt_);
    tree_->SetBranchAddress("genVpt"         , &genVpt_);
    tree_->SetBranchAddress("dmpt"           , &dmpt_);
    tree_->SetBranchAddress("weight"         , &weight_);

    gErrorIgnoreLevel = currentState;
  }

}; 

inline void 
MitLimitTree::InitVariables(){
  // inizialize variables
  mvamet_    = -.1;
  mvametphi_ = -.1;
  jet1pt_    = -.1;
  genjetpt_  = -.1;
  genVpt_    = -.1;
  dmpt_      = -.1;
  weight_    = -.1;

}

#endif

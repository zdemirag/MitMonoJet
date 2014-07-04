#ifndef MitDMSTree_H
#define MitDMSTree_H

#include <set>
#include <vector>
#include <string>
#include <utility>

#include "TFile.h"
#include "TTree.h"
#include "TError.h"
#include "Math/VectorUtil.h"
#include "Math/LorentzVector.h"

class MitDMSTree {
 public:
  /// float doesn't have dictionary by default, so use double
  typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > LorentzVector;

  /// bit map
  /// DON'T CHANGE ORDER
  enum Trigger {
    HLTJetMet   = 1UL<<0,    // event passes jet+met trigger
    HLTMet      = 1UL<<1,    // event passes met trigger
    HLTMuon     = 1UL<<2,    // event passes single muon trigger
    HLTPhoton   = 1UL<<3     // event passes single photon trigger
  };

  /// bit map
  /// DON'T CHANGE ORDER
  enum HLTMatch {
    JetMatch    = 1UL<<0,    // hardest jet is matched to HLT jet
    MuonMatch   = 1UL<<1,    // hardest lepton is matched to HLT muon
    PhotonMatch = 1UL<<2      // hardest photon is matched to HLT photon
  };

  /// bit map
  /// DON'T CHANGE ORDER
  enum Presel {
    Top      = 1UL<<0,    // event passes top preselection
    Wlep     = 1UL<<1,    // event passes W>lv preselection
    Zlep     = 1UL<<2,    // event passes Z>ll preselection
    Met      = 1UL<<3,    // event passes MET preselection
    Vbf      = 1UL<<4,    // event passes VBF preselection
    Gjet     = 1UL<<5     // event passes G+jets preselection
  };

  /// variables
  unsigned int   isData_;
  unsigned int   event_;
  unsigned int   run_;
  unsigned int   lumi_;
  unsigned int   trigger_;
  unsigned int   HLTmatch_;
  unsigned int   nvtx_;
  float          metRaw_;
  float          metRawPhi_;
  float          met_;
  float          metPhi_;
  float          mvamet_;
  float          mvametPhi_;
  float          mvaCov00_;
  float          mvaCov10_;
  float          mvaCov01_;
  float          mvaCov11_;
 
  unsigned int   nlep_;
  LorentzVector  lep1_;
  int            lid1_;
  LorentzVector  lep2_;
  int            lid2_;

  int            ntaus_;
  LorentzVector  tau1_;

  unsigned int   nphotons_;
  LorentzVector  pho1_;

  // here come the tag jets (can contain substructure)
  // substructure jets are only saved for the hardest fat jet
  int            nfjets_;
  LorentzVector  fjet1_;
  float          fjet1CHF_;  
  float          fjet1NHF_;  
  float          fjet1NEMF_; 
  float          fjet1Btag_;
  float          fjet1QGtag_;
  float          fjet1Tau1_;
  float          fjet1Tau2_;
  float          fjet1Tau3_;
  float          fjet1C2b0_;
  float          fjet1C2b0p2_;      
  float          fjet1C2b0p5_;      
  float          fjet1C2b1_;        
  float          fjet1C2b2_;        
  float          fjet1QJetVol_;     
  float          fjet1MassSDb0_;    
  float          fjet1MassSDb2_;    
  float          fjet1MassSDbm1_;   
  float          fjet1MassPruned_;  
  float          fjet1MassFiltered_;
  float          fjet1MassTrimmed_; 
  float          fjet1Pull_;
  float          fjet1PullAngle_;   
  float          fjet1QGtagSub1_;
  float          fjet1QGtagSub2_;
  unsigned int   fjet1PartonId_;
  LorentzVector  fjet2_;
  float          fjet2CHF_;  
  float          fjet2NHF_;  
  float          fjet2NEMF_; 
  float          fjet2Btag_;
  float          fjet2QGtag_;
  float          fjet2Tau1_;
  float          fjet2Tau2_;
  float          fjet2Tau3_;
  float          fjet2C2b0_;
  float          fjet2C2b0p2_;      
  float          fjet2C2b0p5_;      
  float          fjet2C2b1_;        
  float          fjet2C2b2_;        
  float          fjet2QJetVol_;     
  float          fjet2MassSDb0_;    
  float          fjet2MassSDb2_;    
  float          fjet2MassSDbm1_;   
  float          fjet2MassPruned_;  
  float          fjet2MassFiltered_;
  float          fjet2MassTrimmed_; 
  float          fjet2Pull_;
  float          fjet2PullAngle_;   
  float          fjet2QGtagSub1_;
  float          fjet2QGtagSub2_;
  unsigned int   fjet2PartonId_;
 
  unsigned int   nsjets_;
  LorentzVector  sjet1_;
  LorentzVector  sjet2_;

  unsigned int   njets_;
  LorentzVector  jet1_;
  LorentzVector  jet2_;
  LorentzVector  jet3_;
  LorentzVector  jet4_;
  LorentzVector  jet5_;
 
  unsigned int   nbjets_;
  LorentzVector  bjet1_;
  float          bjet1Btag_;
  LorentzVector  bjet2_;
  float          bjet2Btag_;
 
  LorentzVector  genV_;
  unsigned int   genVid_;
  unsigned int   genVdaughterId_;
  float          topPt_;
  float          topBarPt_;
 
  float          Q_;
  float          id1_;
  float          x1_;
  float          pdf1_;
  float          id2_;
  float          x2_;
  float          pdf2_;
  int            processId_;
  float          puweight_;
  float          npu_;
  float          npuPlusOne_;
  float          npuMinusOne_;
  unsigned int   metFiltersWord_;
  unsigned int   preselWord_;

 public:
  /// this is the main element
  TFile *f_;
  TTree *tree_;

  /// hold the names of variables to facilitate things (filled during Init)
  std::vector<std::string> variables_;

  /// default constructor  
  MitDMSTree():
    lepPtr1_(&lep1_),lepPtr2_(&lep2_),
    tauPtr1_(&tau1_),phoPtr1_(&pho1_),
    fjet1Ptr_(&fjet1_),fjet2Ptr_(&fjet2_),
    sjetPtr1_(&sjet1_),sjetPtr2_(&sjet2_),
    jetPtr1_(&jet1_),jetPtr2_(&jet2_),jetPtr3_(&jet3_),jetPtr4_(&jet4_),jetPtr5_(&jet5_),
    bjetPtr1_(&bjet1_),bjetPtr2_(&bjet2_),
    genVPtr_(&genV_) {}
  /// default destructor
  ~MitDMSTree(){
    if (f_) f_->Close();  
  };

  /// initialize varibles and fill list of available variables
  void InitVariables();

  /// load a MitDMSTree
  void LoadTree(const char* file, int type = -1){
    // to load three different ntuples in the same job DMTree0/1
    // type == 0/1 if all variables was added
    // type = -1 (default) if a minimum set of variables was added with tree as name
    f_ = TFile::Open(file);
    assert(f_);
    if     (type == 0) tree_ = dynamic_cast<TTree*>(f_->FindObjectAny("DMSTree"));
    else               tree_ = dynamic_cast<TTree*>(f_->FindObjectAny("tree"));
    assert(tree_);
  }

  /// load a MitDMSTree
  void LoadTree(const char* file, const char* treeName){
    f_ = TFile::Open(file);
    assert(f_);
    tree_ = dynamic_cast<TTree*>(f_->FindObjectAny(treeName));
    assert(tree_);
  }

  /// create a MitDMSTree
  void CreateTree(int type = -1){
    assert(type==type); // just to suppress warnings
    // to create three different ntuples in the same job DMTree0/1
    // type == 0/1 add all variables
    // type = -1 (default) add a minimum set of variables with tree as name
    if     (type == 0) tree_ = new TTree("DMSTree","DM&S ntuples");
    else               tree_ = new TTree("tree","Smurf ntuples");
    f_ = 0;
    InitVariables();
    //book the branches
    tree_->Branch("isData"          , &isData_          ,   "isData/i");
    tree_->Branch("event"           , &event_           ,   "event/i");
    tree_->Branch("run"             , &run_             ,   "run/i");
    tree_->Branch("lumi"            , &lumi_            ,   "lumi/i");
    tree_->Branch("trigger"         , &trigger_         ,   "trigger/i");
    tree_->Branch("HLTmatch"        , &HLTmatch_        ,   "HLTmatch/i");
    tree_->Branch("nvtx"            , &nvtx_            ,   "nvtx/i");
    tree_->Branch("metRaw"          , &metRaw_          ,   "metRaw/F");
    tree_->Branch("metRawPhi"       , &metRawPhi_       ,   "metRawPhi/F");
    tree_->Branch("met"             , &met_             ,   "met/F");
    tree_->Branch("metPhi"          , &metPhi_          ,   "metPhi/F");
    tree_->Branch("mvamet"          , &mvamet_          ,   "mvamet/F");
    tree_->Branch("mvametPhi"       , &mvametPhi_       ,   "mvametPhi/F");
    tree_->Branch("mvaCov00"        , &mvaCov00_        ,   "mvaCov00/F");
    tree_->Branch("mvaCov10"        , &mvaCov10_        ,   "mvaCov10/F");
    tree_->Branch("mvaCov01"        , &mvaCov01_        ,   "mvaCov01/F");
    tree_->Branch("mvaCov11"        , &mvaCov11_        ,   "mvaCov11/F");

    tree_->Branch("nlep"         , &nlep_         ,   "nlep/i");
    tree_->Branch("lep1"         , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >", &lepPtr1_);
    tree_->Branch("lid1"         , &lid1_         ,   "lid1/I");
    tree_->Branch("lep2"         , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >", &lepPtr2_);
    tree_->Branch("lid2"         , &lid2_         ,   "lid2/I");

    tree_->Branch("ntaus"        , &ntaus_        ,   "ntaus/i");
    tree_->Branch("tau1"         , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >", &tauPtr1_);

    tree_->Branch("nphotons"     , &nphotons_     ,   "nphotons/i");
    tree_->Branch("pho1"         , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >", &phoPtr1_);

    tree_->Branch("nfjets", &nfjets_, "nfjets/i");
    tree_->Branch("fjet1", "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >", &fjet1Ptr_);
    tree_->Branch("fjet1CHF"         , &fjet1CHF_         , "fjet1CHF/F");  
    tree_->Branch("fjet1NHF"         , &fjet1NHF_         , "fjet1NHF/F");  
    tree_->Branch("fjet1NEMF"        , &fjet1NEMF_        , "fjet1NEMF/F"); 
    tree_->Branch("fjet1Btag"        , &fjet1Btag_        , "fjet1Btag/F");
    tree_->Branch("fjet1QGtag"       , &fjet1QGtag_       , "fjet1QGtag/F");
    tree_->Branch("fjet1Tau1"        , &fjet1Tau1_        , "fjet1Tau1/F");
    tree_->Branch("fjet1Tau2"        , &fjet1Tau2_        , "fjet1Tau2/F");
    tree_->Branch("fjet1Tau3"        , &fjet1Tau3_        , "fjet1Tau3/F");
    tree_->Branch("fjet1C2b0"        , &fjet1C2b0_        , "fjet1C2b0/F");
    tree_->Branch("fjet1C2b0p2"      , &fjet1C2b0p2_      , "fjet1C2b0p2/F");     
    tree_->Branch("fjet1C2b0p5"      , &fjet1C2b0p5_      , "fjet1C2b0p5/F");     
    tree_->Branch("fjet1C2b1"        , &fjet1C2b1_        , "fjet1C2b1/F");       
    tree_->Branch("fjet1C2b2"        , &fjet1C2b2_        , "fjet1C2b2/F");       
    tree_->Branch("fjet1QJetVol"     , &fjet1QJetVol_     , "fjet1QJetVol/F");    
    tree_->Branch("fjet1MassSDb0"    , &fjet1MassSDb0_    , "fjet1MassSDb0/F");   
    tree_->Branch("fjet1MassSDb2"    , &fjet1MassSDb2_    , "fjet1MassSDb2/F");   
    tree_->Branch("fjet1MassSDbm1"   , &fjet1MassSDbm1_   , "fjet1MassSDbm1/F");  
    tree_->Branch("fjet1MassPruned"  , &fjet1MassPruned_  , "fjet1MassPruned/F"); 
    tree_->Branch("fjet1MassFiltered", &fjet1MassFiltered_, "fjet1MassFiltered/F");
    tree_->Branch("fjet1MassTrimmed" , &fjet1MassTrimmed_ , "fjet1MassTrimmed/F");
    tree_->Branch("fjet1Pull"        , &fjet1Pull_        , "fjet1Pull/F");  
    tree_->Branch("fjet1PullAngle"   , &fjet1PullAngle_   , "fjet1PullAngle/F"); 
    tree_->Branch("fjet1QGtagSub1"   , &fjet1QGtagSub1_   , "fjet1QGtagSub1/F");
    tree_->Branch("fjet1QGtagSub2"   , &fjet1QGtagSub2_   , "fjet1QGtagSub2/F");
    tree_->Branch("fjet1PartonId"    , &fjet1PartonId_    , "fjet1PartonId/i");
    tree_->Branch("fjet2", "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >", &fjet2Ptr_);
    tree_->Branch("fjet2CHF"         , &fjet2CHF_         , "fjet2CHF/F");  
    tree_->Branch("fjet2NHF"         , &fjet2NHF_         , "fjet2NHF/F");  
    tree_->Branch("fjet2NEMF"        , &fjet2NEMF_        , "fjet2NEMF/F"); 
    tree_->Branch("fjet2Btag"        , &fjet2Btag_        , "fjet2Btag/F");
    tree_->Branch("fjet2QGtag"       , &fjet2QGtag_       , "fjet2QGtag/F");
    tree_->Branch("fjet2Tau1"        , &fjet2Tau1_        , "fjet2Tau1/F");
    tree_->Branch("fjet2Tau2"        , &fjet2Tau2_        , "fjet2Tau2/F");
    tree_->Branch("fjet2Tau3"        , &fjet2Tau3_        , "fjet2Tau3/F");
    tree_->Branch("fjet2C2b0"        , &fjet2C2b0_        , "fjet2C2b0/F");
    tree_->Branch("fjet2C2b0p2"      , &fjet2C2b0p2_      , "fjet2C2b0p2/F");     
    tree_->Branch("fjet2C2b0p5"      , &fjet2C2b0p5_      , "fjet2C2b0p5/F");     
    tree_->Branch("fjet2C2b1"        , &fjet2C2b1_        , "fjet2C2b1/F");       
    tree_->Branch("fjet2C2b2"        , &fjet2C2b2_        , "fjet2C2b2/F");       
    tree_->Branch("fjet2QJetVol"     , &fjet2QJetVol_     , "fjet2QJetVol/F");    
    tree_->Branch("fjet2MassSDb0"    , &fjet2MassSDb0_    , "fjet2MassSDb0/F");   
    tree_->Branch("fjet2MassSDb2"    , &fjet2MassSDb2_    , "fjet2MassSDb2/F");   
    tree_->Branch("fjet2MassSDbm1"   , &fjet2MassSDbm1_   , "fjet2MassSDbm1/F");  
    tree_->Branch("fjet2MassPruned"  , &fjet2MassPruned_  , "fjet2MassPruned/F"); 
    tree_->Branch("fjet2MassFiltered", &fjet2MassFiltered_, "fjet2MassFiltered/F");
    tree_->Branch("fjet2MassTrimmed" , &fjet2MassTrimmed_ , "fjet2MassTrimmed/F");
    tree_->Branch("fjet2Pull"        , &fjet2Pull_        , "fjet2Pull/F");  
    tree_->Branch("fjet2PullAngle"   , &fjet2PullAngle_   , "fjet2PullAngle/F"); 
    tree_->Branch("fjet2QGtagSub1"   , &fjet2QGtagSub1_   , "fjet2QGtagSub1/F");
    tree_->Branch("fjet2QGtagSub2"   , &fjet2QGtagSub2_   , "fjet2QGtagSub2/F");
    tree_->Branch("fjet2PartonId"    , &fjet2PartonId_    , "fjet2PartonId/i");
    tree_->Branch("nsjets", &nsjets_, "nsjets/i");
    tree_->Branch("sjet1", "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >", &sjetPtr1_);
    tree_->Branch("sjet2", "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >", &sjetPtr2_);

    tree_->Branch("njets", &njets_, "njets/i");
    tree_->Branch("jet1", "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >", &jetPtr1_);
    tree_->Branch("jet2", "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >", &jetPtr2_);
    tree_->Branch("jet3", "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >", &jetPtr3_);
    tree_->Branch("jet4", "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >", &jetPtr4_);
    tree_->Branch("jet5", "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >", &jetPtr5_);

    tree_->Branch("nbjets", &nbjets_, "nbjets/i");
    tree_->Branch("bjet1", "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >", &bjetPtr1_);
    tree_->Branch("bjet1Btag", &bjet1Btag_, "bjet1Btag/F");
    tree_->Branch("bjet2", "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >", &bjetPtr2_);
    tree_->Branch("bjet2Btag", &bjet2Btag_, "bjet2Btag/F");

    tree_->Branch("genV", "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >", &genVPtr_);
    tree_->Branch("genVid",            &genVid_,             "genVid/i");
    tree_->Branch("genVdaughterId",    &genVdaughterId_,     "genVdaughterId/i");
    tree_->Branch("topPt",             &topPt_    ,          "topPt/F");
    tree_->Branch("topBarPt",          &topBarPt_ ,          "topBarPt/F");

    tree_->Branch("Q",              &Q_	  ,     "Q/F");
    tree_->Branch("id1",            &id1_  ,     "id1/F");
    tree_->Branch("x1",             &x1_	  ,     "x1/F");
    tree_->Branch("pdf1",           &pdf1_ ,     "pdf1/F");
    tree_->Branch("id2",            &id2_  ,     "id2/F");
    tree_->Branch("x2",             &x2_	  ,     "x2/F");
    tree_->Branch("pdf2",           &pdf2_ ,     "pdf2/F");
    tree_->Branch("processId",      &processId_  , "processId/I");
    tree_->Branch("puweight",       &puweight_   , "puweight/F");
    tree_->Branch("npu",            &npu_        , "npu/F");
    tree_->Branch("npuPlusOne",     &npuPlusOne_ , "npuPlusOne/F");
    tree_->Branch("npuMinusOne",    &npuMinusOne_, "npuMinusOne/F");
    tree_->Branch("metFiltersWord", &metFiltersWord_, "metFiltersWord/I");
    tree_->Branch("preselWord",     &preselWord_, "preselWord/I");

  }

  // initialze a MitDMSTree
  void InitTree(int type = -1){
    assert(type==type); // just to suppress warnings
    assert(tree_);
    // don't forget to set pointers to zero before you set address
    // or you will fully appreciate that "ROOT sucks" :)
    InitVariables();
    //Set branch address
    Int_t currentState = gErrorIgnoreLevel;
    // gErrorIgnoreLevel = kError;
    gErrorIgnoreLevel = kBreak;
    tree_->SetBranchAddress("isData",        &isData_);
    tree_->SetBranchAddress("event",         &event_);
    tree_->SetBranchAddress("run",           &run_);
    tree_->SetBranchAddress("lumi",          &lumi_);
    tree_->SetBranchAddress("trigger",       &trigger_);
    tree_->SetBranchAddress("HLTmatch",      &HLTmatch_);
    tree_->SetBranchAddress("nvtx",          &nvtx_);
    tree_->SetBranchAddress("metRaw",        &metRaw_);
    tree_->SetBranchAddress("metRawPhi",     &metRawPhi_);
    tree_->SetBranchAddress("met",           &met_);
    tree_->SetBranchAddress("metPhi",        &metPhi_);
    tree_->SetBranchAddress("mvamet",        &mvamet_);
    tree_->SetBranchAddress("mvametPhi",     &mvametPhi_);
    tree_->SetBranchAddress("mvaCov00",      &mvaCov00_);
    tree_->SetBranchAddress("mvaCov10",      &mvaCov10_);
    tree_->SetBranchAddress("mvaCov01",      &mvaCov01_);
    tree_->SetBranchAddress("mvaCov11",      &mvaCov11_);

    tree_->SetBranchAddress("nlep",          &nlep_);
    tree_->SetBranchAddress("lep1",          &lepPtr1_);
    tree_->SetBranchAddress("lid1",          &lid1_);
    tree_->SetBranchAddress("lep2",          &lepPtr2_);
    tree_->SetBranchAddress("lid2",          &lid2_);

    tree_->SetBranchAddress("ntaus",         &ntaus_);
    tree_->SetBranchAddress("tau1",          &tauPtr1_);

    tree_->SetBranchAddress("nphotons",      &nphotons_);
    tree_->SetBranchAddress("pho1",          &phoPtr1_);

    tree_->SetBranchAddress("nfjets"           , &nfjets_           );
    tree_->SetBranchAddress("fjet1"            , &fjet1Ptr_         );
    tree_->SetBranchAddress("fjet1CHF"         , &fjet1CHF_         );  
    tree_->SetBranchAddress("fjet1NHF"         , &fjet1NHF_         );  
    tree_->SetBranchAddress("fjet1NEMF"        , &fjet1NEMF_        ); 
    tree_->SetBranchAddress("fjet1Btag"        , &fjet1Btag_        );
    tree_->SetBranchAddress("fjet1QGtag"       , &fjet1QGtag_       );
    tree_->SetBranchAddress("fjet1Tau1"        , &fjet1Tau1_        );
    tree_->SetBranchAddress("fjet1Tau2"        , &fjet1Tau2_        );
    tree_->SetBranchAddress("fjet1Tau3"        , &fjet1Tau3_        );
    tree_->SetBranchAddress("fjet1C2b0"        , &fjet1C2b0_        );
    tree_->SetBranchAddress("fjet1C2b0p2"      , &fjet1C2b0p2_      );     
    tree_->SetBranchAddress("fjet1C2b0p5"      , &fjet1C2b0p5_      );     
    tree_->SetBranchAddress("fjet1C2b1"        , &fjet1C2b1_        );       
    tree_->SetBranchAddress("fjet1C2b2"        , &fjet1C2b2_        );       
    tree_->SetBranchAddress("fjet1QJetVol"     , &fjet1QJetVol_     );    
    tree_->SetBranchAddress("fjet1MassSDb0"    , &fjet1MassSDb0_    );   
    tree_->SetBranchAddress("fjet1MassSDb2"    , &fjet1MassSDb2_    );   
    tree_->SetBranchAddress("fjet1MassSDbm1"   , &fjet1MassSDbm1_   );  
    tree_->SetBranchAddress("fjet1MassPruned"  , &fjet1MassPruned_  ); 
    tree_->SetBranchAddress("fjet1MassFiltered", &fjet1MassFiltered_);
    tree_->SetBranchAddress("fjet1MassTrimmed" , &fjet1MassTrimmed_ );
    tree_->SetBranchAddress("fjet1Pull"        , &fjet1Pull_        );  
    tree_->SetBranchAddress("fjet1PullAngle"   , &fjet1PullAngle_   ); 
    tree_->SetBranchAddress("fjet1QGtagSub1"   , &fjet1QGtagSub1_   );
    tree_->SetBranchAddress("fjet1QGtagSub2"   , &fjet1QGtagSub2_   );
    tree_->SetBranchAddress("fjet1PartonId"    , &fjet1PartonId_    );
    tree_->SetBranchAddress("fjet2"            , &fjet2Ptr_         );
    tree_->SetBranchAddress("fjet2CHF"         , &fjet2CHF_         );  
    tree_->SetBranchAddress("fjet2NHF"         , &fjet2NHF_         );  
    tree_->SetBranchAddress("fjet2NEMF"        , &fjet2NEMF_        ); 
    tree_->SetBranchAddress("fjet2Btag"        , &fjet2Btag_        );
    tree_->SetBranchAddress("fjet2QGtag"       , &fjet2QGtag_       );
    tree_->SetBranchAddress("fjet2Tau1"        , &fjet2Tau1_        );
    tree_->SetBranchAddress("fjet2Tau2"        , &fjet2Tau2_        );
    tree_->SetBranchAddress("fjet2Tau3"        , &fjet2Tau3_        );
    tree_->SetBranchAddress("fjet2C2b0"        , &fjet2C2b0_        );
    tree_->SetBranchAddress("fjet2C2b0p2"      , &fjet2C2b0p2_      );     
    tree_->SetBranchAddress("fjet2C2b0p5"      , &fjet2C2b0p5_      );     
    tree_->SetBranchAddress("fjet2C2b1"        , &fjet2C2b1_        );       
    tree_->SetBranchAddress("fjet2C2b2"        , &fjet2C2b2_        );       
    tree_->SetBranchAddress("fjet2QJetVol"     , &fjet2QJetVol_     );    
    tree_->SetBranchAddress("fjet2MassSDb0"    , &fjet2MassSDb0_    );   
    tree_->SetBranchAddress("fjet2MassSDb2"    , &fjet2MassSDb2_    );   
    tree_->SetBranchAddress("fjet2MassSDbm1"   , &fjet2MassSDbm1_   );  
    tree_->SetBranchAddress("fjet2MassPruned"  , &fjet2MassPruned_  ); 
    tree_->SetBranchAddress("fjet2MassFiltered", &fjet2MassFiltered_);
    tree_->SetBranchAddress("fjet2MassTrimmed" , &fjet2MassTrimmed_ );
    tree_->SetBranchAddress("fjet2Pull"        , &fjet2Pull_        );  
    tree_->SetBranchAddress("fjet2PullAngle"   , &fjet2PullAngle_   ); 
    tree_->SetBranchAddress("fjet2QGtagSub1"   , &fjet2QGtagSub1_   );
    tree_->SetBranchAddress("fjet2QGtagSub2"   , &fjet2QGtagSub2_   );
    tree_->SetBranchAddress("fjet2PartonId"    , &fjet2PartonId_    );

    tree_->SetBranchAddress("nsjets"           , &nsjets_           );
    tree_->SetBranchAddress("sjet1"            , &sjetPtr1_         );
    tree_->SetBranchAddress("sjet2"            , &sjetPtr2_         );

    tree_->SetBranchAddress("njets"            , &njets_            );
    tree_->SetBranchAddress("jet1"             , &jetPtr1_          );
    tree_->SetBranchAddress("jet2"             , &jetPtr2_          );
    tree_->SetBranchAddress("jet3"             , &jetPtr3_          );
    tree_->SetBranchAddress("jet4"             , &jetPtr4_          );
    tree_->SetBranchAddress("jet5"             , &jetPtr5_          );

    tree_->SetBranchAddress("nbjets"           , &nbjets_           );
    tree_->SetBranchAddress("bjet1"            , &bjetPtr1_         );
    tree_->SetBranchAddress("bjet1Btag"        , &bjet1Btag_        );
    tree_->SetBranchAddress("bjet2"            , &bjetPtr2_         );
    tree_->SetBranchAddress("bjet2Btag"        , &bjet2Btag_        );

    tree_->SetBranchAddress("genV"             , &genVPtr_          );
    tree_->SetBranchAddress("genVid "          , &genVid_           );
    tree_->SetBranchAddress("genVdaughterId"   , &genVdaughterId_   );
    tree_->SetBranchAddress("topPt"            , &topPt_            );
    tree_->SetBranchAddress("topBarPt"         , &topBarPt_         );

    tree_->SetBranchAddress("Q"             ,	&Q_             );
    tree_->SetBranchAddress("id1"           ,	&id1_           );
    tree_->SetBranchAddress("x1"            ,	&x1_            );
    tree_->SetBranchAddress("pdf1"          ,	&pdf1_          );
    tree_->SetBranchAddress("id2"           ,	&id2_           );
    tree_->SetBranchAddress("x2"            ,	&x2_            );
    tree_->SetBranchAddress("pdf2"          ,	&pdf2_          );
    tree_->SetBranchAddress("processId"     , &processId_     );
    tree_->SetBranchAddress("puweight"      , &puweight_      );
    tree_->SetBranchAddress("npu"           ,	&npu_           );
    tree_->SetBranchAddress("npuPlusOne"    , &npuPlusOne_    );
    tree_->SetBranchAddress("npuMinusOne"   , &npuMinusOne_   );
    tree_->SetBranchAddress("metFiltersWord" , &metFiltersWord_ );
    tree_->SetBranchAddress("preselWord"    , &preselWord_);

    gErrorIgnoreLevel = currentState;
  }

  private:

  LorentzVector* lepPtr1_;
  LorentzVector* lepPtr2_;
  LorentzVector* tauPtr1_;
  LorentzVector* phoPtr1_;

  LorentzVector* fjet1Ptr_;
  LorentzVector* fjet2Ptr_;
  LorentzVector* sjetPtr1_;
  LorentzVector* sjetPtr2_;

  LorentzVector* jetPtr1_;
  LorentzVector* jetPtr2_;
  LorentzVector* jetPtr3_;
  LorentzVector* jetPtr4_;
  LorentzVector* jetPtr5_;

  LorentzVector* bjetPtr1_;
  LorentzVector* bjetPtr2_;

  LorentzVector* genVPtr_;
}; 

inline void 
MitDMSTree::InitVariables(){
  // inizialize variables
  isData_        = 0;
  event_         = 0;
  run_           = 0;
  lumi_          = 0;
  trigger_       = 0;
  HLTmatch_      = 0;
  nvtx_          = 0;
  metRaw_        = -999.;
  metRawPhi_     = -999.;
  met_           = -999.;
  metPhi_        = -999.;
  mvamet_        = -999.;
  mvametPhi_     = -999.;
  mvaCov00_      = -999.;
  mvaCov10_      = -999.; 
  mvaCov01_      = -999.;
  mvaCov11_      = -999.;

  nlep_          = 0;
  lep1_       	 = LorentzVector();
  lid1_          = 0;
  lep2_       	 = LorentzVector();
  lid2_          = 0;

  ntaus_         = 0;
  tau1_       	 = LorentzVector();
  
  nphotons_      = 0;
  pho1_       	 = LorentzVector();

  nfjets_         = 0;
  fjet1_          = LorentzVector();
  fjet1CHF_       = -999.;
  fjet1NHF_       = -999.;
  fjet1NEMF_      = -999.;
  fjet1Btag_      = -999.;
  fjet1QGtag_     = -999.;
  fjet1Tau1_      = -999.;
  fjet1Tau2_      = -999.;
  fjet1Tau3_      = -999.;
  fjet1C2b0_      = -999.;
  fjet1C2b0p2_      = -999.;
  fjet1C2b0p5_      = -999.;
  fjet1C2b1_        = -999.;
  fjet1C2b2_        = -999.;
  fjet1QJetVol_     = -999.;
  fjet1MassSDb0_    = -999.;
  fjet1MassSDb2_    = -999.;
  fjet1MassSDbm1_   = -999.;
  fjet1MassPruned_  = -999.;
  fjet1MassFiltered_= -999.;
  fjet1MassTrimmed_ = -999.;
  fjet1Pull_        = -999.;
  fjet1PullAngle_   = -999.;
  fjet1QGtagSub1_   = -999.;
  fjet1QGtagSub2_   = -999.;
  fjet1PartonId_  = 0;
  fjet2_          = LorentzVector();
  fjet2CHF_       = -999.;
  fjet2NHF_       = -999.;
  fjet2NEMF_      = -999.;
  fjet2Btag_      = -999.;
  fjet2QGtag_     = -999.;
  fjet2Tau1_      = -999.;
  fjet2Tau2_      = -999.;
  fjet2Tau3_      = -999.;
  fjet2C2b0_      = -999.;
  fjet2C2b0p2_      = -999.;
  fjet2C2b0p5_      = -999.;
  fjet2C2b1_        = -999.;
  fjet2C2b2_        = -999.;
  fjet2QJetVol_     = -999.;
  fjet2MassSDb0_    = -999.;
  fjet2MassSDb2_    = -999.;
  fjet2MassSDbm1_   = -999.;
  fjet2MassPruned_  = -999.;
  fjet2MassFiltered_= -999.;
  fjet2MassTrimmed_ = -999.;
  fjet2Pull_        = -999.;
  fjet2PullAngle_   = -999.;
  fjet2QGtagSub1_   = -999.;
  fjet2QGtagSub2_   = -999.;
  fjet2PartonId_  = 0;
  
  nsjets_        = 0;
  sjet1_         = LorentzVector();
  sjet2_         = LorentzVector();
  
  njets_         = 0;
  jet1_          = LorentzVector();
  jet2_          = LorentzVector();
  jet3_          = LorentzVector();
  jet4_          = LorentzVector();
  jet5_          = LorentzVector();

  nbjets_        = 0; 
  bjet1_         = LorentzVector();
  bjet1Btag_     = 0;
  bjet2_         = LorentzVector();
  bjet2Btag_     = 0;

  genV_          = LorentzVector();
  genVid_        = 0;
  genVdaughterId_= 0;
  topPt_         = 0;
  topBarPt_      = 0;
  
  Q_		         = -999.;
  id1_  	       = -999.;
  x1_		         = -999.;
  pdf1_ 	       = -999.;  
  id2_  	       = -999.;  
  x2_		         = -999.;
  pdf2_ 	       = -999.;  
  processId_	   = 0;
  puweight_      = 1.;
  npu_           = -999.;
  npuPlusOne_    = -999.;
  npuMinusOne_   = -999.;
  metFiltersWord_= 0;  
  preselWord_    = 0;  
}

#endif

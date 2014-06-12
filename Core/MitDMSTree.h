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

  // here comes the tag jet (can contain substructure)
  LorentzVector  tjet_;
  float          tjetCHF_;  
  float          tjetNHF_;  
  float          tjetNEMF_; 
  float          tjetBtag_;
  float          tjetQGtag_;
  float          tjetTau1_;
  float          tjetTau2_;
  float          tjetTau3_;
  float          tjetC2b0_;
  float          tjetC2b0p2_;      
  float          tjetC2b0p5_;      
  float          tjetC2b1_;        
  float          tjetC2b2_;        
  float          tjetQJetVol_;     
  float          tjetMassSDb0_;    
  float          tjetMassSDb2_;    
  float          tjetMassSDbm1_;   
  float          tjetMassPruned_;  
  float          tjetMassFiltered_;
  float          tjetMassTrimmed_; 
  unsigned int   tjetPartonId_;
 
  unsigned int   nsjets_;
  LorentzVector  sjet1_;
  LorentzVector  sjet2_;

  unsigned int   njets_;
  LorentzVector  jet1_;
  float          jet1Btag_;
  LorentzVector  jet2_;
  float          jet2Btag_;
  LorentzVector  jet3_;
  float          jet3Btag_;
  LorentzVector  jet4_;
  float          jet4Btag_;
  LorentzVector  jet5_;
  float          jet5Btag_;
  unsigned int   nbjets_;
 
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
    tjetPtr_(&tjet_),
    sjetPtr1_(&sjet1_),sjetPtr2_(&sjet2_),
    jetPtr1_(&jet1_),jetPtr2_(&jet2_),jetPtr3_(&jet3_),jetPtr4_(&jet4_),jetPtr5_(&jet5_){}
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
    if     (type == 0) tree_ = dynamic_cast<TTree*>(f_->Get("DMSTree"));
    else	             tree_ = dynamic_cast<TTree*>(f_->Get("tree"));
    assert(tree_);
  }

  /// load a MitDMSTree
  void LoadTree(const char* file, const char* treeName){
    f_ = TFile::Open(file);
    assert(f_);
    tree_ = dynamic_cast<TTree*>(f_->Get(treeName));
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

    tree_->Branch("tjet"          , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >", &tjetPtr_);
    tree_->Branch("tjetCHF"       , &tjetCHF_      ,   "tjetCHF/F");  
    tree_->Branch("tjetNHF"       , &tjetNHF_      ,   "tjetNHF/F");  
    tree_->Branch("tjetNEMF"      , &tjetNEMF_     ,   "tjetNEMF/F"); 
    tree_->Branch("tjetBtag"      , &tjetBtag_     ,   "tjetBtag/F");
    tree_->Branch("tjetQGtag"     , &tjetQGtag_    ,   "tjetQGtag/F");
    tree_->Branch("tjetTau1"      , &tjetTau1_     ,   "tjetTau1/F");
    tree_->Branch("tjetTau2"      , &tjetTau2_     ,   "tjetTau2/F");
    tree_->Branch("tjetTau3"      , &tjetTau3_     ,   "tjetTau3/F");
    tree_->Branch("tjetC2b0"      , &tjetC2b0_     ,   "tjetC2b0/F");
    tree_->Branch("tjetC2b0p2"       , &tjetC2b0p2_         , "tjetC2b0p2/F");     
    tree_->Branch("tjetC2b0p5"       , &tjetC2b0p5_         , "tjetC2b0p5/F");     
    tree_->Branch("tjetC2b1"         , &tjetC2b1_           , "tjetC2b1/F");       
    tree_->Branch("tjetC2b2"         , &tjetC2b2_           , "tjetC2b2/F");       
    tree_->Branch("tjetQJetVol"      , &tjetQJetVol_        , "tjetQJetVol/F");    
    tree_->Branch("tjetMassSDb0"     , &tjetMassSDb0_       , "tjetMassSDb0/F");   
    tree_->Branch("tjetMassSDb2"     , &tjetMassSDb2_       , "tjetMassSDb2/F");   
    tree_->Branch("tjetMassSDbm1"    , &tjetMassSDbm1_      , "tjetMassSDbm1/F");  
    tree_->Branch("tjetMassPruned"   , &tjetMassPruned_     , "tjetMassPruned/F"); 
    tree_->Branch("tjetMassFiltered" , &tjetMassFiltered_   , "tjetMassFiltered/F");
    tree_->Branch("tjetMassTrimmed"  , &tjetMassTrimmed_    , "tjetMassTrimmed/F");
    tree_->Branch("tjetPartonId"  , &tjetPartonId_ ,   "tjetPartonId/i");

    tree_->Branch("nsjets"        , &nsjets_        ,   "nsjets/i");
    tree_->Branch("sjet1"         , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >", &sjetPtr1_);
    tree_->Branch("sjet2"         , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >", &sjetPtr2_);

    tree_->Branch("njets"         , &njets_         ,   "njets/i");
    tree_->Branch("jet1"          , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >", &jetPtr1_);
    tree_->Branch("jet1Btag"      , &jet1Btag_      ,   "jet1Btag/F");
    tree_->Branch("jet2"          , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >", &jetPtr2_);
    tree_->Branch("jet2Btag"      , &jet2Btag_      ,   "jet2Btag/F");
    tree_->Branch("jet3"          , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >", &jetPtr3_);
    tree_->Branch("jet3Btag"      , &jet3Btag_      ,   "jet3Btag/F");
    tree_->Branch("jet4"          , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >", &jetPtr4_);
    tree_->Branch("jet4Btag"      , &jet4Btag_      ,   "jet4Btag/F");
    tree_->Branch("jet5"          , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >", &jetPtr5_);
    tree_->Branch("jet5Btag"      , &jet5Btag_      ,   "jet5Btag/F");
    tree_->Branch("nbjets"        , &nbjets_        ,   "nbjets/i");

    tree_->Branch("Q",             &Q_	  ,     "Q/F");
    tree_->Branch("id1",           &id1_  ,     "id1/F");
    tree_->Branch("x1",            &x1_	  ,     "x1/F");
    tree_->Branch("pdf1",          &pdf1_ ,     "pdf1/F");
    tree_->Branch("id2",           &id2_  ,     "id2/F");
    tree_->Branch("x2",            &x2_	  ,     "x2/F");
    tree_->Branch("pdf2",          &pdf2_ ,     "pdf2/F");
    tree_->Branch("processId",     &processId_  , "processId/I");
    tree_->Branch("puweight",      &puweight_   , "puweight/F");
    tree_->Branch("npu",           &npu_        , "npu/F");
    tree_->Branch("npuPlusOne",    &npuPlusOne_ , "npuPlusOne/F");
    tree_->Branch("npuMinusOne",   &npuMinusOne_, "npuMinusOne/F");
    tree_->Branch("metFiltersWord",&metFiltersWord_, "metFiltersWord/I");
    tree_->Branch("preselWord",    &preselWord_, "preselWord/I");

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

    tree_->SetBranchAddress("ntaus"        , &ntaus_);
    tree_->SetBranchAddress("tau1"         , &tauPtr1_);

    tree_->SetBranchAddress("nphotons"                  , &nphotons_);
    tree_->SetBranchAddress("pho1"                      , &phoPtr1_);

    tree_->SetBranchAddress("tjet"          , &tjetPtr_       );
    tree_->SetBranchAddress("tjetCHF"       , &tjetCHF_       );  
    tree_->SetBranchAddress("tjetNHF"       , &tjetNHF_       );  
    tree_->SetBranchAddress("tjetNEMF"      , &tjetNEMF_      ); 
    tree_->SetBranchAddress("tjetBtag"      , &tjetBtag_      );
    tree_->SetBranchAddress("tjetQGtag"     , &tjetQGtag_     );
    tree_->SetBranchAddress("tjetTau1"      , &tjetTau1_      );
    tree_->SetBranchAddress("tjetTau2"      , &tjetTau2_      );
    tree_->SetBranchAddress("tjetTau3"      , &tjetTau3_      );
    tree_->SetBranchAddress("tjetC2b0"      , &tjetC2b0_      );
    tree_->SetBranchAddress("tjetC2b0p2"       , &tjetC2b0p2_      );
    tree_->SetBranchAddress("tjetC2b0p5"       , &tjetC2b0p5_      );
    tree_->SetBranchAddress("tjetC2b1"         , &tjetC2b1_        );
    tree_->SetBranchAddress("tjetC2b2"         , &tjetC2b2_        );
    tree_->SetBranchAddress("tjetQJetVol"      , &tjetQJetVol_     );
    tree_->SetBranchAddress("tjetMassSDb0"     , &tjetMassSDb0_    );
    tree_->SetBranchAddress("tjetMassSDb2"     , &tjetMassSDb2_    );
    tree_->SetBranchAddress("tjetMassSDbm1"    , &tjetMassSDbm1_   );
    tree_->SetBranchAddress("tjetMassPruned"   , &tjetMassPruned_  );
    tree_->SetBranchAddress("tjetMassFiltered" , &tjetMassFiltered_);
    tree_->SetBranchAddress("tjetMassTrimmed"  , &tjetMassTrimmed_ );
    tree_->SetBranchAddress("tjetPartonId"  , &tjetPartonId_  );

    tree_->SetBranchAddress("nsjets"        , &nsjets_        );
    tree_->SetBranchAddress("sjet1"         , &sjetPtr1_      );
    tree_->SetBranchAddress("sjet2"         , &sjetPtr2_      );

    tree_->SetBranchAddress("njets"         , &njets_         );
    tree_->SetBranchAddress("jet1"          , &jetPtr1_       );
    tree_->SetBranchAddress("jet1Btag"      , &jet1Btag_      );
    tree_->SetBranchAddress("jet2"          , &jetPtr2_       );
    tree_->SetBranchAddress("jet2Btag"      , &jet2Btag_      );
    tree_->SetBranchAddress("jet3"          , &jetPtr3_       );
    tree_->SetBranchAddress("jet3Btag"      , &jet3Btag_      );
    tree_->SetBranchAddress("jet4"          , &jetPtr4_       );
    tree_->SetBranchAddress("jet4Btag"      , &jet4Btag_      );
    tree_->SetBranchAddress("jet5"          , &jetPtr5_       );
    tree_->SetBranchAddress("jet5Btag"      , &jet5Btag_      );
    tree_->SetBranchAddress("nbjets"        , &nbjets_        );

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

  LorentzVector* tjetPtr_;
  LorentzVector* sjetPtr1_;
  LorentzVector* sjetPtr2_;
  LorentzVector* jetPtr1_;
  LorentzVector* jetPtr2_;
  LorentzVector* jetPtr3_;
  LorentzVector* jetPtr4_;
  LorentzVector* jetPtr5_;
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

  tjet_          = LorentzVector();
  tjetCHF_       = -999.;
  tjetNHF_       = -999.;
  tjetNEMF_      = -999.;
  tjetBtag_      = -999.;
  tjetQGtag_     = -999.;
  tjetTau1_      = -999.;
  tjetTau2_      = -999.;
  tjetTau3_      = -999.;
  tjetC2b0_      = -999.;
  tjetC2b0p2_      = -999.;
  tjetC2b0p5_      = -999.;
  tjetC2b1_        = -999.;
  tjetC2b2_        = -999.;
  tjetQJetVol_     = -999.;
  tjetMassSDb0_    = -999.;
  tjetMassSDb2_    = -999.;
  tjetMassSDbm1_   = -999.;
  tjetMassPruned_  = -999.;
  tjetMassFiltered_= -999.;
  tjetMassTrimmed_ = -999.;
  tjetPartonId_  = 0;
  
  nsjets_        = 0;
  sjet1_         = LorentzVector();
  sjet2_         = LorentzVector();
  
  njets_         = 0;
  jet1_          = LorentzVector();
  jet1Btag_      = 0;
  jet2_          = LorentzVector();
  jet2Btag_      = 0;
  jet3_          = LorentzVector();
  jet3Btag_      = 0; 
  jet4_          = LorentzVector();
  jet4Btag_      = 0; 
  jet5_          = LorentzVector();
  jet5Btag_      = 0; 
  nbjets_        = 0; 
  
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

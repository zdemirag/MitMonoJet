 // $Id $

#include "MitMonoJet/SelMods/interface/MonoJetAnalysisMod.h"
#include <TH1D.h>
#include <TH2D.h>
#include <TParameter.h>
#include "MitAna/DataUtil/interface/Debug.h"
#include "MitAna/DataTree/interface/Names.h"
#include "MitPhysics/Init/interface/ModNames.h"
#include "MitAna/DataCont/interface/ObjArray.h"
#include "MitCommon/MathTools/interface/MathUtils.h"
#include "MitPhysics/Utils/interface/JetTools.h"
#include "MitPhysics/Utils/interface/MetTools.h"
#include "MitAna/DataTree/interface/ParticleCol.h"
#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include <iostream>
#include <sstream>

using namespace mithep;
ClassImp(mithep::MonoJetAnalysisMod)

//--------------------------------------------------------------------------------------------------
MonoJetAnalysisMod::MonoJetAnalysisMod(const char *name, const char *title) : 
  BaseMod(name,title),
  // define all the Branches to load
  fMetBranchName                 ("PFMet"),
  fJetsName			             (Names::gkPFJetBrn), //added by TJ
  fMetFromBranch                 (kTRUE),
  fJetsFromBranch                (kTRUE),
  // ----------------------------------------
  // collections....
  fMet                           (0),
  fJets				             (0), //added by TJ
  fMinNumJets		 	         (1),
  fMinJetEt			             (30),
  fMaxJetEta			         (2.4),
  fMinMetEt 			         (30),  
  // counters....
  fNEventsSelected               (0)
{
  // Constructor.
}

//--------------------------------------------------------------------------------------------------
void MonoJetAnalysisMod::Begin()
{
  // Run startup code on the client machine. For this module, we dont do
  // anything here.
}

//--------------------------------------------------------------------------------------------------
void MonoJetAnalysisMod::SlaveBegin()
{
  // Run startup code on the computer (slave) doing the actual analysis. Here,
  // we typically initialize histograms and other analysis objects and request
  // branches. For this module, we request a branch of the MitTree.

  // Load Branches
  ReqEventObject(fMetBranchName,   fMet,      fMetFromBranch);
  ReqEventObject(fJetsName,        fJets,	  fJetsFromBranch); 

  //Create your histograms here
  //*************************************************************************************************
  // Selection Histograms
  //*************************************************************************************************
  AddTH1(fMonoJetSelection,"hMonoJetSelection", ";Cut Number;Number of Events",             17, -1.5, 15.5);

  // Create const char* labels with cut values
  std::ostringstream os1;
  std::string s1;
  os1<<"N Jets >= "<<fMinNumJets;
  s1 = os1.str();
  const char* labelCut1 = s1.c_str();
  std::ostringstream os2;
  std::string s2;
  os2<<"Jet Et >= "<<fMinJetEt;
  s2 = os2.str();
  const char* labelCut2 = s2.c_str();
  std::ostringstream os3;
  std::string s3;
  os3<<"Jet Eta <= "<<fMaxJetEta;
  s3 = os3.str();
  const char* labelCut3 = s3.c_str();
  std::ostringstream os4;
  std::string s4;
  os4<<"Met >= "<<fMinMetEt;
  s4 = os4.str();
  const char* labelCut4 = s4.c_str();

  // Set selection histogram bin labels
  fMonoJetSelection->GetXaxis()->TAxis::SetBinLabel(1, "All Events");
  fMonoJetSelection->GetXaxis()->TAxis::SetBinLabel(2, labelCut1);
  fMonoJetSelection->GetXaxis()->TAxis::SetBinLabel(3, labelCut2);
  fMonoJetSelection->GetXaxis()->TAxis::SetBinLabel(4, labelCut3);
  fMonoJetSelection->GetXaxis()->TAxis::SetBinLabel(5, labelCut4);

  //***********************************************************************************************
  // Histograms after preselection
  //***********************************************************************************************
  AddTH1(fMetEt              ,"hMetEt", ";MetEt;Number of Events" ,400,0.,400.);
  AddTH1(fJetEt		         ,"hJetEt", ";JetEt;Number of Events" ,400,0.,400.); //added by TJ
  AddTH1(fJetEta	         ,"hJetEta",";JetEta;Number of Events",50,-5,5); //added by TJ for testing

}

//--------------------------------------------------------------------------------------------------
void MonoJetAnalysisMod::Process()
{
  // Process entries of the tree.
  LoadEventObject(fMetBranchName,   fMet);
  assert(fMet);
  LoadEventObject(fJetsName,   fJets); //added by TJ
  assert(fJets);  
  //cout << "I am able to load event object!" << "\n";


  //*********************************************************************************************
  //Define Cuts (I changed this from 3 to 4 to include eta - TJ)
  //*********************************************************************************************
  const int nCuts = 4;
  bool passCut[nCuts] = {
	  false, 
	  false, 
	  false,
	  false,
	  };
	  
  //***********************************************************************************************
  //Discard events with no identified photons (jets? - TJ)
  //***********************************************************************************************
  // if (fPhotons->GetEntries() > 0)  passCut[0] = true;
  if (fJets->GetEntries() > fMinNumJets) passCut[0] = true;

  //***********************************************************************************************
  //Discard events with soft jets - (changed by TJ)
  //***********************************************************************************************
  //added by TJ
  int nHardJet = 0;
  for (UInt_t i = 0; i < fJets->GetEntries(); ++i) {
	  const Jet *jet = fJets->At(i);
	  if (jet->Et() <= fMinJetEt ) continue;
	  nHardJet++;
  }
  if (nHardJet > 0) passCut[1] = true;

  //**********************************************************************************************
  //Discard events that are outside of eta range (added by TJ)
  //*********************************************************************************************
  int nHardEtaJet = 0;
  for (UInt_t i = 0; i < fJets->GetEntries(); ++i) {
	const Jet *jet = fJets->At(i);
	if (jet->Et() <= fMinJetEt ) continue;
	if (TMath::Abs(jet->Eta()) >= fMaxJetEta) continue;
	nHardEtaJet++;
  }
  if (nHardEtaJet > 0) passCut[2] = true;
  
  

  //***********************************************************************************************
  //Discard events with soft met
  //***********************************************************************************************
  const Met *stdMet = fMet->At(0);
  if ( stdMet->Pt() > fMinMetEt )  passCut[3] = true;
	  
  //*********************************************************************************************
  //Make Selection Histograms. Number of events passing each level of cut
  //*********************************************************************************************  
  //Cut Selection Histograms
  fMonoJetSelection->Fill(-1,1);

  for (int k=0;k<nCuts;k++) {
    bool pass = true;
    bool passPreviousCut = true;
    for (int p=0;p<=k;p++) {
      pass = (pass && passCut[p]);
      if (p<k)
        passPreviousCut = (passPreviousCut&& passCut[p]);
    }
    if (pass)
      fMonoJetSelection->Fill(k,1);
  }

  
  bool passAllCuts = true;
  for(int c=0; c<nCuts; c++) passAllCuts = passAllCuts & passCut[c];
  if(passAllCuts) {
	  fNEventsSelected++;
      //*****************************************************************************************
	  //Make Preselection Histograms  
	  //*****************************************************************************************
	  fMetEt->Fill(fMet->At(0)->Et());
	  fJetEt->Fill(fJets->At(0)->Et());
	  fJetEta->Fill(fJets->At(0)->Eta());
  }
  else 
	  this->SkipEvent(); //skip the event if does not passes the cuts
  
  return;
}

//--------------------------------------------------------------------------------------------------
void MonoJetAnalysisMod::SlaveTerminate()
{
  
  // Run finishing code on the computer (slave) that did the analysis. For this
  // module, we dont do anything here.
  cout << "selected events on MonoJetAnalysisMod: " << fNEventsSelected << endl;

} 
//--------------------------------------------------------------------------------------------------
void MonoJetAnalysisMod::Terminate()
{
  // Run finishing code on the client computer. For this module, we dont do
  // anything here.

}

// $Id: HLTEvtSelMod.cc,v 1.3 2010/11/12 10:31:40 ceballos Exp $

#include "MitMonoJet/Mods/interface/HLTEvtSelMod.h"
#include <TH1D.h>
#include <TH2D.h>
#include <TParameter.h>
#include <TTree.h>
#include "MitPhysics/Init/interface/ModNames.h"
#include "MitAna/DataCont/interface/ObjArray.h"
#include "MitAna/DataTree/interface/PFMet.h"
#include "MitCommon/MathTools/interface/MathUtils.h"
#include "MitAna/DataTree/interface/CompositeParticleCol.h"
#include "MitAna/DataTree/interface/MetCol.h"
#include "MitAna/DataTree/interface/GenericParticleCol.h"
#include "MitAna/DataTree/interface/JetCol.h"
#include "MitAna/DataTree/interface/VertexCol.h"
#include "TDataMember.h"
#include "TFile.h"

using namespace mithep;
ClassImp(mithep::HLTEvtSelMod)

//--------------------------------------------------------------------------------------------------
HLTEvtSelMod::HLTEvtSelMod(const char *name, const char *title) : 
  BaseMod(name,title),
  fMetName("PFMet"),
  fEvtHdrName(Names::gkEvtHeaderBrn),
  fVertexName(Names::gkPVBeamSpotBrn),
  fVertices(0),
  fEventHeader(0),
  fTupleName("hMonoJetTreeHLT")
{
  // Constructor.
}

//--------------------------------------------------------------------------------------------------
void HLTEvtSelMod::Begin()
{
  // Run startup code on the client machine. For this module, we dont do
  // anything here.  
}

//--------------------------------------------------------------------------------------------------
void HLTEvtSelMod::Process()
{
  // Process entries of the tree. For this module, we just load the branches and  
  //Obtain all the good objects from the event cleaning module
  JetOArr *CleanJets           = GetObjThisEvt<JetOArr>(fCleanJetsName);

  //Get the MET branch
  LoadEventObject(fMetName,           fMet,           true);

  fVertices = GetObjThisEvt<VertexOArr>(fVertexName);
  LoadBranch(fEvtHdrName);
  
  const TriggerObjectCol *objs = GetHLTObjects(fObjsName);

  if (!objs){
    printf("HLTEvtSelMod::TriggerObjectCol not found\n");
    return;
  }

  // Prepare all the bits necessary to get the HLT-offline match information

  bool HLTMu  = false;

  // Jet matching, one entry for each offline jet
  bool HLTJets[4]  = {false,false,false,false};
  // Jet no matching, just store the trigger decision
  bool HLTJetsNoMatch  = false;

  // MET : first bit is for CaloMET, second for PFMHT
  bool HLTMets[2]  = {false,false};

  fMonoJetEventHLT->fMETEt = fMet->At(0)->Pt();
    
  // Loop on trigger obects
  for(unsigned int iHLT=0; iHLT<objs->GetEntries(); ++iHLT) {

    const TriggerObject* to = objs->At(iHLT);
      
    // Muon HLT
    if (to->Type() == TriggerObject::TriggerMuon && TString(to->TrigName()).Contains("IsoMu24_eta2p1_v"))
      HLTMu = true;

    // jet HLT matching, loop on jets and check if they are matched to the relevant HLT jets
    fMonoJetEventHLT->fNjets = CleanJets->GetEntries();
    for (unsigned int iJet = 0; iJet < CleanJets->GetEntries(); iJet++) {
      fMonoJetEventHLT->fJetPt[iJet] = CleanJets->At(iJet)->Mom().Pt();
      fMonoJetEventHLT->fJetEta[iJet] = CleanJets->At(iJet)->Mom().Eta();
      if (to->Type() == TriggerObject::TriggerJet && MathUtils::DeltaR(CleanJets->At(iJet)->Mom(), to->Mom()) < 0.5) {
        if (TString(to->TrigName()).Contains("MonoCentralPFJet80_PFMETnoMu95_NHEF0p95_v")) HLTJets[iJet] = true;
        if (TString(to->TrigName()).Contains("MonoCentralPFJet80_PFMETnoMu105_NHEF0p95_v")) HLTJets[iJet] = true;
      }
    }// end loop on jets  
    if (to->Type() == TriggerObject::TriggerJet) {
      if (TString(to->TrigName()).Contains("MonoCentralPFJet80_PFMETnoMu95_NHEF0p95_v") || TString(to->TrigName()).Contains("MonoCentralPFJet80_PFMETnoMu105_NHEF0p95_v")) HLTJetsNoMatch = true;
    }

    // MET HLT
    if (to->Type() == TriggerObject::TriggerMET && TString(to->TrigName()).Contains("HLT_MET120_HBHENoiseCleaned_v"))
      HLTMets[0] = true;
    if ( to->Type() == TriggerObject::TriggerMHT && 
        (TString(to->TrigName()).Contains("MonoCentralPFJet80_PFMETnoMu95_NHEF0p95_v") || TString(to->TrigName()).Contains("MonoCentralPFJet80_PFMETnoMu95_NHEF0p95_v")) 
       )
      HLTMets[1] = true;


  }// end loop on HLT obj

  // fill the ntuple variables
  fMonoJetEventHLT->fHLTMu = HLTMu;
  fMonoJetEventHLT->fHLTJets[0] = HLTJets[0];
  fMonoJetEventHLT->fHLTJets[1] = HLTJets[1];
  fMonoJetEventHLT->fHLTJets[2] = HLTJets[2];
  fMonoJetEventHLT->fHLTJets[3] = HLTJets[3];
  fMonoJetEventHLT->fHLTJetsNoMatch = HLTJetsNoMatch;
  fMonoJetEventHLT->fHLTMets[0] = HLTMets[0];
  fMonoJetEventHLT->fHLTMets[1] = HLTMets[1];

  fMonoJetEventHLT->fNumVertices = fVertices->GetEntries();
  fMonoJetEventHLT->fRunNumber = fEventHeader->RunNum();

  fMonoJetTupleHLT->Fill();
}
//--------------------------------------------------------------------------------------------------
void HLTEvtSelMod::SlaveBegin()
{
  // Run startup code on the computer (slave) doing the actual analysis. Here,
  // we typically initialize histograms and other analysis objects and request
  // branches. For this module, we request a branch of the MitTree.
  ReqEventObject(fMetName,    fMet,            true);
  ReqBranch(fEvtHdrName,      fEventHeader);

  //***********************************************************************************************
  //Create Ntuple Tree  
  //***********************************************************************************************
  fMonoJetEventHLT = new MonoJetEventHLT;
  TFile *ftmp = TFile::Open(TString::Format("%s_tmp.root",GetName()),"RECREATE");
  fMonoJetTupleHLT = new TTree(fTupleName.Data(),fTupleName.Data());
    
  //make flattish tree from classes so we don't have to rely on dictionaries for reading later
  TClass *eclass = TClass::GetClass("mithep::MonoJetEventHLT");
  TList  *elist  = eclass->GetListOfDataMembers();

  for (int i=0; i<elist->GetEntries(); ++i) {
    const TDataMember *tdm = static_cast<const TDataMember*>(elist->At(i));//ming
    if (!(tdm->IsBasic() && tdm->IsPersistent())) continue;
    TString typestring;
    if (TString(tdm->GetTypeName()).BeginsWith("Char_t")) typestring = "B";
    else if (TString(tdm->GetTypeName()).BeginsWith("UChar_t")) typestring = "b";
    else if (TString(tdm->GetTypeName()).BeginsWith("Short_t")) typestring = "S";
    else if (TString(tdm->GetTypeName()).BeginsWith("UShort_t")) typestring = "s";
    else if (TString(tdm->GetTypeName()).BeginsWith("Int_t")) typestring = "I";
    else if (TString(tdm->GetTypeName()).BeginsWith("UInt_t")) typestring = "i";
    else if (TString(tdm->GetTypeName()).BeginsWith("Float_t")) typestring = "F";
    else if (TString(tdm->GetTypeName()).BeginsWith("Double_t")) typestring = "D";
    else if (TString(tdm->GetTypeName()).BeginsWith("Long64_t")) typestring = "L";
    else if (TString(tdm->GetTypeName()).BeginsWith("ULong64_t")) typestring = "l";
    else if (TString(tdm->GetTypeName()).BeginsWith("Bool_t")) typestring = "O";
    else continue;
    //printf("%s %s: %i\n",tdm->GetTypeName(),tdm->GetName(),int(tdm->GetOffset()));
    Char_t *addr = (Char_t*)fMonoJetEventHLT;//ming:?
    assert(sizeof(Char_t)==1);
    if (TString(tdm->GetName()).CompareTo("fHLTJets") == 0) fMonoJetTupleHLT->Branch(tdm->GetName(),addr + tdm->GetOffset(),TString::Format("%s[4]/%s",tdm->GetName(),typestring.Data())); 
    else if (TString(tdm->GetName()).CompareTo("fHLTMets") == 0) fMonoJetTupleHLT->Branch(tdm->GetName(),addr + tdm->GetOffset(),TString::Format("%s[2]/%s",tdm->GetName(),typestring.Data())); 
    else if (TString(tdm->GetName()).CompareTo("fJetPt") == 0) fMonoJetTupleHLT->Branch(tdm->GetName(),addr + tdm->GetOffset(),TString::Format("%s[4]/%s",tdm->GetName(),typestring.Data())); 
    else if (TString(tdm->GetName()).CompareTo("fJetEta") == 0) fMonoJetTupleHLT->Branch(tdm->GetName(),addr + tdm->GetOffset(),TString::Format("%s[4]/%s",tdm->GetName(),typestring.Data())); 
    else fMonoJetTupleHLT->Branch(tdm->GetName(),addr + tdm->GetOffset(),TString::Format("%s/%s",tdm->GetName(),typestring.Data()));
  }
  
  AddOutput(fMonoJetTupleHLT);
}

//--------------------------------------------------------------------------------------------------
void HLTEvtSelMod::SlaveTerminate()
{
  // Run finishing code on the computer (slave) that did the analysis
}

//--------------------------------------------------------------------------------------------------
void HLTEvtSelMod::Terminate()
{
  // Run finishing code on the client computer
}

//--------------------------------------------------------------------------------------------------
// $Id: HltEvtSelMod.h,v 1.1 2010/10/23 04:44:54 ceballos Exp $
//
// HltEvtSelMod
//
//
// Authors: L. Di Matteo, C. Paus
//--------------------------------------------------------------------------------------------------
#ifndef MITMONOJET_HltEvtSelMod_H
#define MITMONOJET_HltEvtSelMod_H

#include "MitAna/TreeMod/interface/BaseMod.h" 
#include "MitAna/DataTree/interface/CollectionsFwd.h"

namespace mithep 
{
  class MonoJetEventHlt
  {
  public:
    // trigger bits per object
    Bool_t  fHltMu;           // Muon HLT, used to unbias the sample
    Bool_t  fHltJets[4];      // HLT jet matching, one per offline jet ((0-3)
    Bool_t  fHltJetsNoMatch;  // Jet no matching trigger decision
    Bool_t  fHltMets[2];      // MET : first bit is for CaloMET, second for PFMHT
    // offline quantities
    Float_t fMetEt;           // Missing Et
    Int_t   fNjets;           // number of jets
    Float_t fJetPt[4];        // jet pt (0-3)
    Float_t fJetEta[4];       // jet eta (0-3)
    Int_t   fNumVertices;     // number of vertices (primary?)
    Int_t   fRunNumber;       // run number

    void    ResetTriggerBits()
    {
      // Muon trigger
      fHltMu = false;
      // Matched jet trigger
      fHltJets[0] = false; fHltJets[1] = false; fHltJets[2] = false; fHltJets[3] = false;
      // Unmatched jet trigger decision
      fHltJetsNoMatch = false;
      // MET: first bit is for CaloMET, second for PFMHT
      fHltMets[0] = false; fHltMets[1] = false;
    }
  };
  
  class HltEvtSelMod : public BaseMod
  {
  public:
    HltEvtSelMod(const char *name  = "HltEvtSelMod", 
		 const char *title = "Module to select the HLT triggers");
    ~HltEvtSelMod() {}

    // Defining the names of the collections
    void                 SetCleanJetsName(TString s)     { fCleanJetsName = s;}	
    void                 SetTrigObjsName (const char *n) { fHltObjsName = n;  }
    void                 SetMetName      (TString s)     { fMetName = s;      }   
    void                 SetVertexName   (TString s)     { fVertexName = s;   }   
    // Cut variables
    void                 SetPtJetCut     (double x)      { fPtJetCut = x;     }
    void                 SetEtaJetCut    (double x)      { fEtaJetCut = x;    }
    		         
  protected:	         
    // Standard methods of a module
    void                 Begin();
    void                 Process();
    void                 SlaveBegin();
    void                 SlaveTerminate();
    void                 Terminate();      

  private:
    // A little helper
    TString              TypeString(const TDataMember *tdm);

    // Event header
    TString              fEvtHdrName;      // name of event header (this is available by default?!)
    const EventHeader   *fEventHeader;     // event header (do we really need this?)
    // Vertex collextion
    TString              fVertexName;      // name of vertex collection
    const VertexCol     *fVertices;        // vertex collection
    // HL Trigger objects
    TString              fHltObjsName;     // name of the HLT objects to consider
    // MET name
    TString              fMetName;         // name of particle flow MET collection
    const PFMetCol      *fMet;             // particle flow MET collection
    // Cleaned jets
    TString              fCleanJetsName;   // name of the offline jet collection
    // Cuts
    Double_t             fTrgMatchDeltaR;  // deltaR matching between offline and trigger jet
    Double_t             fPtJetCut;        // minimum Pt cut on jet
    Double_t             fEtaJetCut;       // maximum Eta cut on jet
    
    // What we store
    TString              fTupleName;           // name o the ntuple
    MonoJetEventHlt     *fMonoJetEventHlt;     // event trigger summary (relevent for us)
    TTree               *fMonoJetTupleHlt;     // the monojet information
    
    ClassDef(HltEvtSelMod,1) // Analysis class definition
  };
}
#endif

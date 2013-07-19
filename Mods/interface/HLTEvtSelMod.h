//------------------------------------------------------------------------------
// $Id: HLTEvtSelMod.h,v 1.1 2010/10/23 04:44:54 ceballos Exp $
//
// HLTEvtSelMod
//
//
// Authors: ceballos
//------------------------------------------------------------------------------

#ifndef MITMONOJET_HLTEvtSelMod_H
#define MITMONOJET_HLTEvtSelMod_H

#include "MitAna/TreeMod/interface/BaseMod.h" 
#include "MitAna/DataTree/interface/CollectionsFwd.h"
#include "TRandom.h"

class TH1D;
class TH2D;

namespace mithep 
{
  class MonoJetEventHLT
  {
    public:  
      // Muon HLT, used to unbias the sample
      Bool_t fHLTMu;
      // Jet matching, one entry for each offline jet
      Bool_t fHLTJets[4];
      // Jet no matching, just store the trigger decision
      Bool_t fHLTJetsNoMatch;
      // MET : first bit is for CaloMET, second for PFMHT
      Bool_t fHLTMets[2];
      // offline quantities
      Float_t fMETEt;
      Int_t   fNjets;
      Float_t fJetPt[4];
      Float_t fJetEta[4];
      Int_t fNumVertices;
      Int_t fRunNumber;
  };
  
  class HLTEvtSelMod : public BaseMod
  {
    public:
    HLTEvtSelMod(const char *name="HLTEvtSelMod", 
		 const char *title="Example analysis module with all branches");
      ~HLTEvtSelMod() {}
      void         SetPrInt_tDebug(Bool_t b)     { fPrInt_tDebug = b; }	
      void         SetCleanJetsName (TString s)  { fCleanJetsName = s;}	
      void         SetTrigObjsName(const char *n){ fObjsName = n;     }
      void         SetMetName(TString s)         { fMetName = s;      }   
      void         SetVertexName(TString s)      { fVertexName = s;   }   
      void         SetPtJetCut(double x)         { fPtJetCut = x;     }
      void         SetEtaJetCut(double x)        { fEtaJetCut = x;    }

    protected:
      Bool_t      fPrInt_tDebug;
      Double_t    fPtJetCut;
      Double_t    fEtaJetCut;
      TString     fMetName;
      TString     fCleanJetsName;
      TString             fEvtHdrName;
      const EventHeader  *fEventHeader;
      TString             fVertexName;
      const VertexCol    *fVertices;
      TString      fObjsName;
      const PFMetCol     *fMet;

      void         Begin();
      void         Process();
      void         SlaveBegin();
      void         SlaveTerminate();
      void         Terminate();      

      TString           fTupleName;
      MonoJetEventHLT  *fMonoJetEventHLT;
      TTree            *fMonoJetTupleHLT;
  
      ClassDef(HLTEvtSelMod,1) // TAM example analysis module
  };
}
#endif

//--------------------------------------------------------------------------------------------------
// $Id: PFJet.h,v 1.7 2012/03/28 12:15:34 paus Exp $
//
// PFJet
//
// This class holds information about reconstructed jets and their substructure
//
// Authors: L.DiMatteo
//--------------------------------------------------------------------------------------------------

#ifndef MITMONOJET_DATATREE_XLFATJET_H
#define MITMONOJET_DATATREE_XLFATJET_H
 
#include "MitAna/DataTree/interface/Jet.h"
#include "MitAna/DataCont/interface/RefArray.h"
#include "MitMonoJet/DataTree/interface/XlSubJet.h"

namespace mithep 
{
  class XlFatJet : public Jet
  {
    public:
      XlFatJet() : 
                fTau1(0), fTau2(0), fTau3(0) {}
      XlFatJet(Double_t px, Double_t py, Double_t pz, Double_t e) : 
                Jet(px,py,pz,e),
                fTau1(0), fTau2(0), fTau3(0) {}
      XlFatJet(const Jet & p) : 
                Jet(p),
                fTau1(0), fTau2(0), fTau3(0) {}

      void                  AddSubJet(const XlSubJet *p)          { fSubJets.Add(p); ClearCharge();}
      Bool_t                HasSubJet(const XlSubJet *p)    const { return fSubJets.HasObject(p);  }
      Jet                  *MakeCopy()                      const { return new XlFatJet(*this);    }
      UInt_t                NSubJets()                      const { return fSubJets.Entries();     }
      Double_t              Tau1()                          const { return fTau1;                  }
      Double_t              Tau2()                          const { return fTau2;                  }
      Double_t              Tau3()                          const { return fTau3;                  }
      EObjType              ObjType()                       const { return kXlFatJet;              } 
      const XlSubJet       *SubJet(UInt_t i)                const { return fSubJets.At(i);         }
      void                  SetTau1(Double_t t)                   { fTau1     = t;                 }
      void                  SetTau2(Double_t t)                   { fTau2     = t;                 }
      void                  SetTau3(Double_t t)                   { fTau3     = t;                 }

      // Some structural tools
      void                  Mark(UInt_t i=1)                const;

    protected:

      Double32_t            fTau1;      //1-subjettiness
      Double32_t            fTau2;      //2-subjettiness
      Double32_t            fTau3;      //3-subjettiness
      RefArray<XlSubJet>    fSubJets;   //sub jets in the jet

    ClassDef(XlFatJet, 1) // XlFatJet class
  };
}

//--------------------------------------------------------------------------------------------------
inline void mithep::XlFatJet::Mark(UInt_t ib) const
{
  // mark myself
  mithep::DataObject::Mark(ib);
  // mark my dependencies if they are there
  fSubJets.Mark(ib);
}

#endif

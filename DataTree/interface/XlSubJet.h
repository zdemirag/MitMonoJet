//--------------------------------------------------------------------------------------------------
// $Id: PFJet.h,v 1.7 2012/03/28 12:15:34 paus Exp $
//
// PFJet
//
// This class holds information about reconstructed sub-jets
//
// Authors: L.DiMatteo
//--------------------------------------------------------------------------------------------------

#ifndef MITMONOJET_DATATREE_XLSUBJET_H
#define MITMONOJET_DATATREE_XLSUBJET_H
 
#include "MitAna/DataTree/interface/Jet.h"
#include "MitAna/DataTree/interface/Types.h"
#include "MitCommon/DataFormats/interface/Vect3.h"

namespace mithep 
{
  class XlSubJet : public Jet
  {
    public:
      enum ESubJetType {
        eX = 0,      //unidentified
        eV,          //vector boson subjet
        eTop         //top subjet
      };


      XlSubJet() : 
                fBTag(0), fQGTag(0), fSubJetType(eX) {}
      XlSubJet(Double_t px, Double_t py, Double_t pz, Double_t e) : 
                Jet(px,py,pz,e),
                fBTag(0), fQGTag(0), fSubJetType(eX) {}
      XlSubJet(const Jet & p) : 
                Jet(p),
                fBTag(0), fQGTag(0), fSubJetType(eX) {}

      Jet                  *MakeCopy()                      const { return new XlSubJet(*this);    }
      ThreeVector           Axis()                          const { return fAxis.V();              }
      Double_t              BTag()                          const { return fBTag;                  }
      Double_t              QGTag()                         const { return fQGTag;                 }
      ESubJetType           SubJetType()                    const { return fSubJetType;            } 
      void                  SetAxis(const ThreeVector &a)         { fAxis     = a;                 }
      void                  SetBTag(Double_t d)                   { fBTag     = d;                 }
      void                  SetQGTag(Double_t d)                  { fQGTag    = d;                 }
      void                  SetSubJetType(ESubJetType t)          { fSubJetType = t;               } 

      // Some structural tools
      void                  Mark(UInt_t i=1)                const;

    protected:

      Vect3                 fAxis;       //sub-jet axis found by the jet subclustering algorithm
      Double32_t            fBTag;       //sub-jet b-tag value
      Double32_t            fQGTag;      //sub-jet quark-gluon tag value
      ESubJetType           fSubJetType; //subjet type
      
    ClassDef(XlSubJet, 1) // XlFatJet class
  };
}

//--------------------------------------------------------------------------------------------------
inline void mithep::XlSubJet::Mark(UInt_t ib) const
{
  // mark myself
  mithep::DataObject::Mark(ib);
}

#endif

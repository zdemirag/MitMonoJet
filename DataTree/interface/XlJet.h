//--------------------------------------------------------------------------------------------------
// $Id: PFJet.h,v 1.7 2012/03/28 12:15:34 paus Exp $
//
// XlJet
//
// This class holds information about reconstructed jets and some extra info
//
// Authors: L.DiMatteo
//--------------------------------------------------------------------------------------------------

#ifndef MITMONOJET_DATATREE_XLJET_H
#define MITMONOJET_DATATREE_XLJET_H
 
#include "MitAna/DataTree/interface/PFJet.h"
#include "MitAna/DataTree/interface/Types.h"
#include "MitCommon/DataFormats/interface/Vect4M.h"

namespace mithep 
{
  class XlJet : public PFJet
  {
    public:
      XlJet() :
                fCharge (0),
                fQGTag(0), 
                fPullY(0), fPullPhi(0) {}                 
      XlJet(Double_t px, Double_t py, Double_t pz, Double_t e) : 
                PFJet(px,py,pz,e),
                fCharge (0),
                fQGTag(0), 
                fPullY(0), fPullPhi(0) {}                 
      XlJet(const PFJet & p) : 
                PFJet(p),
                fCharge (0),
                fQGTag(0), 
                fPullY(0), fPullPhi(0) {}                 

      Double_t              Charge()                        const { return fCharge;                } 
      Double_t              QGTag()                         const { return fQGTag;                 } 
      Double_t              PullY()                         const { return fPullY;                 }
      Double_t              PullPhi()                       const { return fPullPhi;               }

      void                  SetCharge()                           { fCharge  = this->GetCharge();  } 
      void                  SetQGTag(Double_t t)                  { fQGTag   = t;                  } 
      void                  SetPullY(Double_t t)                  { fPullY   = t;                  }
      void                  SetPullPhi(Double_t t)                { fPullPhi = t;                  }

      // Some structural tools
      void                  Mark(UInt_t i=1)                const;

    protected:
      Double32_t            GetCharge()                     const; 
      
      Double32_t            fCharge;       //Pt-weighted jet charge
      Double32_t            fQGTag;        //QG tagging
      Double32_t            fPullY;        //Color pull Y directon
      Double32_t            fPullPhi;      //Color pull phi directon

    ClassDef(XlJet, 0) // XlJet class
  };
}

//--------------------------------------------------------------------------------------------------
inline void mithep::XlJet::Mark(UInt_t ib) const
{
  // mark myself
  mithep::DataObject::Mark(ib);
}

//--------------------------------------------------------------------------------------------------
inline Double32_t mithep::XlJet::GetCharge() const
{
  // Return the sum of constituents PFCandidates weighted by their pt
  
  Double_t charge = 0;
  for (UInt_t i=0; i<NPFCands(); ++i)
    charge += PFCand(i)->Charge()*PFCand(i)->Pt();

  return charge/this->Pt();
}

#endif

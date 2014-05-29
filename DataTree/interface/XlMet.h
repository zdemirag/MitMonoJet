//--------------------------------------------------------------------------------------------------
// $Id: PF.h,v 1.7 2012/03/28 12:15:34 paus Exp $
//
// PF
//
// This class holds information about MVA met
//
// Authors: L.DiMatteo
//--------------------------------------------------------------------------------------------------

#ifndef MITMONOJET_DATATREE_XLMET_H
#define MITMONOJET_DATATREE_XLMET_H
 
#include "MitAna/DataTree/interface/Met.h"
#include <TMatrixD.h>

namespace mithep 
{
  class XlMet : public Met
  {
    public:

      XlMet() : 
                fMetCov(2,2) {}
      XlMet(Double_t mex, Double_t mey) : 
                Met(mex, mey), fMetCov(2,2) {}

      Met                  *MakeCopy()                      const { return new XlMet(*this);       }
      TMatrixD              GetCovMatrix()                  const { return fMetCov;                }
      double                Cov00()                         const { return fMetCov(0,0);           }
      double                Cov10()                         const { return fMetCov(1,0);           }
      double                Cov01()                         const { return fMetCov(0,1);           }
      double                Cov11()                         const { return fMetCov(1,1);           }
      void                  SetCovMatrix(const TMatrixD &m)       { fMetCov   = m;                 }

      // Some structural tools
      void                  Mark(UInt_t i=1)                const;

    protected:
      TMatrixD              fMetCov;     //MET covariance Matrix
      
    ClassDef(XlMet, 1) // XlMet class
  };
}

//--------------------------------------------------------------------------------------------------
inline void mithep::XlMet::Mark(UInt_t ib) const
{
  // mark myself
  mithep::DataObject::Mark(ib);
}

#endif

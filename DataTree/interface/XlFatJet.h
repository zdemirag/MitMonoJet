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
 
#include "MitAna/DataTree/interface/PFJet.h"
#include "MitAna/DataCont/interface/RefArray.h"
#include "MitAna/DataTree/interface/Types.h"
#include "MitCommon/DataFormats/interface/Vect4M.h"
#include "MitMonoJet/DataTree/interface/XlSubJet.h"

namespace mithep 
{
  class XlFatJet : public PFJet
  {
    public:
      XlFatJet() :
                fQGTag(0), 
                fTau1(0), fTau2(0), fTau3(0),
                fC2b0(0), fC2b0p2(0), fC2b0p5(0), fC2b1(0), fC2b2(0),
                fQJetVol(0), 
                fMassSDb0(0), fMassSDb2(0), fMassSDbm1(0),
                fMassPruned(0), fMassFiltered(0), fMassTrimmed(0) {}                 
      XlFatJet(Double_t px, Double_t py, Double_t pz, Double_t e) : 
                PFJet(px,py,pz,e),
                fQGTag(0), 
                fTau1(0), fTau2(0), fTau3(0),
                fC2b0(0), fC2b0p2(0), fC2b0p5(0), fC2b1(0), fC2b2(0),
                fQJetVol(0), 
                fMassSDb0(0), fMassSDb2(0), fMassSDbm1(0),
                fMassPruned(0), fMassFiltered(0), fMassTrimmed(0) {}                 
      XlFatJet(const PFJet & p) : 
                PFJet(p),
                fQGTag(0), 
                fTau1(0), fTau2(0), fTau3(0),
                fC2b0(0), fC2b0p2(0), fC2b0p5(0), fC2b1(0), fC2b2(0),
                fQJetVol(0), 
                fMassSDb0(0), fMassSDb2(0), fMassSDbm1(0),
                fMassPruned(0), fMassFiltered(0), fMassTrimmed(0) {}                 

      void                  AddSubJet(const XlSubJet *p)          { fSubJets.Add(p);               }
      Bool_t                HasSubJet(const XlSubJet *p)    const { return fSubJets.HasObject(p);  }
      Jet                  *MakeCopy()                      const { return new XlFatJet(*this);    }
      UInt_t                NSubJets()                      const { return fSubJets.Entries();     }
      Double_t              QGTag()                         const { return fQGTag;                 } 
      Double_t              Tau1()                          const { return fTau1;                  }
      Double_t              Tau2()                          const { return fTau2;                  }
      Double_t              Tau3()                          const { return fTau3;                  }
      Double_t              C2b0()                          const { return fC2b0;                  }
      Double_t              C2b0p2()                        const { return fC2b0p2;                }
      Double_t              C2b0p5()                        const { return fC2b0p5;                }
      Double_t              C2b1()                          const { return fC2b1;                  }
      Double_t              C2b2()                          const { return fC2b2;                  }
      Double_t              QJetVol()                       const { return fQJetVol;               }
      Double_t              MassSDb0()                      const { return fMassSDb0;              }
      Double_t              MassSDb2()                      const { return fMassSDb2;              }
      Double_t              MassSDbm1()                     const { return fMassSDbm1;             }
      Double_t              MassPruned()                    const { return fMassPruned;            }
      Double_t              MassFiltered()                  const { return fMassFiltered;          }
      Double_t              MassTrimmed()                   const { return fMassTrimmed;           }

      const XlSubJet       *SubJet(UInt_t i)                const { return fSubJets.At(i);         }
      void                  SetQGTag(Double_t t)                  { fQGTag       = t;              } 
      void                  SetTau1(Double_t t)                   { fTau1        = t;              }
      void                  SetTau2(Double_t t)                   { fTau2        = t;              }
      void                  SetTau3(Double_t t)                   { fTau3        = t;              }
      void                  SetC2b0(Double_t t)                   { fC2b0        = t;              }
      void                  SetC2b0p2(Double_t t)                 { fC2b0p2      = t;              }
      void                  SetC2b0p5(Double_t t)                 { fC2b0p5      = t;              }
      void                  SetC2b1(Double_t t)                   { fC2b1        = t;              }
      void                  SetC2b2(Double_t t)                   { fC2b2        = t;              }
      void                  SetQJetVol(Double_t t)                { fQJetVol     = t;              }
      void                  SetMassSDb0(Double_t t)               { fMassSDb0    = t;              }
      void                  SetMassSDb2(Double_t t)               { fMassSDb2    = t;              }
      void                  SetMassSDbm1(Double_t t)              { fMassSDbm1   = t;              }
      void                  SetMassPruned(Double_t t)             { fMassPruned  = t;              }
      void                  SetMassFiltered(Double_t t)           { fMassFiltered = t;             }
      void                  SetMassTrimmed(Double_t t)            { fMassTrimmed = t;              }

      // Some structural tools
      void                  Mark(UInt_t i=1)                const;

    protected:

      Double32_t            fQGTag;        //QG tagging
      Double32_t            fTau1;         //1-subjettiness
      Double32_t            fTau2;         //2-subjettiness
      Double32_t            fTau3;         //3-subjettiness
      Double32_t            fC2b0;         //ECF ratio order 2, beta 0
      Double32_t            fC2b0p2;       //ECF ratio order 2, beta 0.2
      Double32_t            fC2b0p5;       //ECF ratio order 2, beta 0.2
      Double32_t            fC2b1;         //ECF ratio order 2, beta 1
      Double32_t            fC2b2;         //ECF ratio order 2, beta 2      
      Double32_t            fQJetVol;      //QJets volatility
      Double32_t            fMassSDb0;     //Groomed mass (soft drop b 0)
      Double32_t            fMassSDb2;     //Groomed mass (soft drop b 2)
      Double32_t            fMassSDbm1;    //Groomed mass (soft drop b-1)
      Double32_t            fMassPruned;   //Groomed mass (pruning)
      Double32_t            fMassFiltered; //Groomed mass (filtering)
      Double32_t            fMassTrimmed;  //Groomed mass (trimming)
      RefArray<XlSubJet>    fSubJets;      //sub jets in the jet

    ClassDef(XlFatJet, 2) // XlFatJet class
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



//--------------------------------------------------------------------------------------------------
// $Id: PFJet.h,v 1.7 2012/03/28 12:15:34 paus Exp $
//
// XsIsoParticle
//
// This class holds information about reconstructed leptons and their quality criteria parameters 
// (isolation,id)
//
// Authors: L.DiMatteo
//--------------------------------------------------------------------------------------------------

#ifndef MITMONOJET_DATATREE_XSISOPARTICLE_H
#define MITMONOJET_DATATREE_XSISOPARTICLE_H
 
#include "MitAna/DataTree/interface/Particle.h"
#include "MitCommon/DataFormats/interface/Vect4M.h"

namespace mithep 
{
  class XsIsoParticle : public Particle
  {
    public:
      enum EParticleId {
        eX = 0,      //unidentified
        eTightMuon   //passing tight muon identification
      };

      XsIsoParticle() :
                fCharge(0),
                fChHadIso(0),fNeHadIso(0),fGammaIso(0),
                fParticleId(eX) {}
      XsIsoParticle(Double_t px, Double_t py, Double_t pz, Double_t e) : 
                fMom(FourVector(px,py,pz,e)),
                fCharge(0),
                fChHadIso(0),fNeHadIso(0),fGammaIso(0),
                fParticleId(eX) {}
      XsIsoParticle(const Particle & p) : 
                Particle(p),
                fCharge(0),
                fChHadIso(0),fNeHadIso(0),fGammaIso(0),
                fParticleId(eX) {}

      Double_t              ChHadIso()                      const { return fChHadIso;          }
      Double_t              NeHadIso()                      const { return fNeHadIso;          }
      Double_t              GammaIso()                      const { return fGammaIso;          }
      EParticleId           ParticleId()                    const { return fParticleId;        } 

      void                  SetCharge(Double_t t)                 { fCharge = t; ClearCharge();}
      void                  SetChHadIso(Double_t t)               { fChHadIso        = t;      }
      void                  SetNeHadIso(Double_t t)               { fNeHadIso        = t;      }
      void                  SetGammaIso(Double_t t)               { fGammaIso        = t;      }
      void                  SetParticleId(EParticleId t)          { fParticleId      = t;      } 
      void                  SetMom(Double_t px, Double_t py, Double_t pz, Double_t e);
 
      // Some structural tools
      void                  Mark(UInt_t i=1)                const;

    protected:
      Double_t              GetCharge()                     const;
      void                  GetMom()                        const;
   
      Vect4M                fMom;          //four momentum vector      
      Double32_t            fCharge;       //particle em charge
      Double32_t            fChHadIso;     //PF Charged Hadron isolation
      Double32_t            fNeHadIso;     //PF Neutral Hadron isolation
      Double32_t            fGammaIso;     //PF Photon isolation
      EParticleId           fParticleId;   //Iso particle id flag

    ClassDef(XsIsoParticle, 0) // XsIsoParticle class
  };
}

//--------------------------------------------------------------------------------------------------
inline void mithep::XsIsoParticle::Mark(UInt_t ib) const
{
  // mark myself
  mithep::DataObject::Mark(ib);
}

//--------------------------------------------------------------------------------------------------
inline Double_t mithep::XsIsoParticle::GetCharge() const
{
  // Get stored charge
  return fCharge;
}

//--------------------------------------------------------------------------------------------------
inline void mithep::XsIsoParticle::GetMom() const
{
  // Get momentum values from stored values.
  fCachedMom.SetCoordinates(fMom.Pt(),fMom.Eta(),fMom.Phi(),fMom.M()); 
}

//--------------------------------------------------------------------------------------------------
inline void mithep::XsIsoParticle::SetMom(Double_t px, Double_t py, Double_t pz, Double_t e)
{ 
  // Set four vector.
  fMom.SetXYZT(px, py, pz, e);
  ClearMom();
}

#endif

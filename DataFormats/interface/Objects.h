#ifndef monojet_Objects_h
#define monojet_Objects_h
#include "Rtypes.h"

namespace monojet {

  class JetCollection;
  class MetCollection;
  class ElectronCollection;
  class MuonCollection;
  class TauCollection;
  class PhotonCollection;

  class Jet {
  public:
    Jet(JetCollection&, UInt_t idx);

    Float_t& pt;
    Float_t& eta;
    Float_t& phi;
    Float_t& btag;
    Float_t& chFrac;
    Float_t& nhFrac;
    Float_t& neFrac;
  };

  class Met {
  public:
    Met(MetCollection&, UInt_t idx);

    Float_t& met;
    Float_t& phi;
  };

  class Electron {
  public:
    Electron(ElectronCollection&, UInt_t idx);

    Float_t& pt;
    Float_t& eta;
    Float_t& phi;
  };

  class Muon {
  public:
    Muon(MuonCollection&, UInt_t idx);

    Float_t& pt;
    Float_t& eta;
    Float_t& phi;
  };

  class Tau {
  public:
    Tau(TauCollection&, UInt_t idx);

    Float_t& pt;
    Float_t& eta;
    Float_t& phi;
    Float_t& disc;
  };

  class Photon {
  public:
    Photon(PhotonCollection&, UInt_t idx);

    Float_t& pt;
    Float_t& eta;
    Float_t& phi;
  };

}

#endif

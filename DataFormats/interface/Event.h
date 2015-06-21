#ifndef monojet_Event_h
#define monojet_Event_h
#include "MitMonoJet/DataFormats/interface/Collections.h"

namespace monojet {

  class Event {
  public:
    UInt_t runNumber;
    UInt_t lumiNumber;
    UInt_t eventNumber;

    JetCollection jets;
    MetCollection mets;
    ElectronCollection electrons;
    MuonCollection muons;
    TauCollection taus;
    PhotonCollection photons;

    void setAddress(TTree&);
    void book(TTree&);
  };

}

#endif

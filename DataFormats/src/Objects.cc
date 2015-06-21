#include "MitMonoJet/DataFormats/interface/Objects.h"

#include "MitMonoJet/DataFormats/interface/Collections.h"

monojet::Jet::Jet(JetCollection& col, UInt_t idx) :
  pt(col.pt[idx]),
  eta(col.eta[idx]),
  phi(col.phi[idx]),
  btag(col.btag[idx]),
  chFrac(col.chFrac[idx]),
  nhFrac(col.nhFrac[idx]),
  neFrac(col.neFrac[idx])
{
}

monojet::Met::Met(MetCollection& col, UInt_t idx) :
  met(col.met[idx]),
  phi(col.phi[idx])
{
}

monojet::Electron::Electron(ElectronCollection& col, UInt_t idx) :
  pt(col.pt[idx]),
  eta(col.eta[idx]),
  phi(col.phi[idx])
{
}

monojet::Muon::Muon(MuonCollection& col, UInt_t idx) :
  pt(col.pt[idx]),
  eta(col.eta[idx]),
  phi(col.phi[idx])
{
}

monojet::Tau::Tau(TauCollection& col, UInt_t idx) :
  pt(col.pt[idx]),
  eta(col.eta[idx]),
  phi(col.phi[idx]),
  disc(col.disc[idx])
{
}

monojet::Photon::Photon(PhotonCollection& col, UInt_t idx) :
  pt(col.pt[idx]),
  eta(col.eta[idx]),
  phi(col.phi[idx])
{
}


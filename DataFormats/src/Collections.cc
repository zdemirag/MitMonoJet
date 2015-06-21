#include "MitMonoJet/DataFormats/interface/Collections.h"
#include "TTree.h"
#include <stdexcept>
#include <memory>

monojet::JetCollection::JetCollection() :
  array_(std::allocator<Jet>().allocate(NMAX))
{
  for (unsigned iP(0); iP != NMAX; ++iP)
    new (array_ + iP) Jet(*this, iP);
}

monojet::JetCollection::~JetCollection()
{
  std::allocator<Jet>().deallocate(array_, NMAX);
}

monojet::JetCollection::reference
monojet::JetCollection::at(UInt_t _idx)
{
  if (_idx < size)
    return array_[_idx];
  throw std::out_of_range("JetCollection::at");
}

monojet::JetCollection::const_reference
monojet::JetCollection::at(UInt_t _idx) const
{
  if (_idx < size)
    return array_[_idx];
  throw std::out_of_range("JetCollection::at");
}

void
monojet::JetCollection::resize(UInt_t _size)
{
  if (size < NMAX) {
    size = _size;
    return;
  }
  throw std::length_error("JetCollection::resize");
}

void
monojet::JetCollection::setAddress(TTree& _tree)
{
  _tree.SetBranchAddress("jet.size", &size);
  _tree.SetBranchAddress("jet.pt", pt);
  _tree.SetBranchAddress("jet.eta", eta);
  _tree.SetBranchAddress("jet.phi", phi);
  _tree.SetBranchAddress("jet.btag", btag);
  _tree.SetBranchAddress("jet.chFrac", chFrac);
  _tree.SetBranchAddress("jet.nhFrac", nhFrac);
  _tree.SetBranchAddress("jet.neFrac", neFrac);
}

void
monojet::JetCollection::book(TTree& _tree)
{
  _tree.Branch("jet.size", &size, "size/i");
  _tree.Branch("jet.pt", pt, "pt[jet.size]/F");
  _tree.Branch("jet.eta", eta, "eta[jet.size]/F");
  _tree.Branch("jet.phi", phi, "phi[jet.size]/F");
  _tree.Branch("jet.btag", btag, "btag[jet.size]/F");
  _tree.Branch("jet.chFrac", chFrac, "chFrac[jet.size]/F");
  _tree.Branch("jet.nhFrac", nhFrac, "nhFrac[jet.size]/F");
  _tree.Branch("jet.neFrac", neFrac, "neFrac[jet.size]/F");
}

monojet::MetCollection::MetCollection() :
  array_(std::allocator<Met>().allocate(NMAX))
{
  for (unsigned iP(0); iP != NMAX; ++iP)
    new (array_ + iP) Met(*this, iP);
}

monojet::MetCollection::~MetCollection()
{
  std::allocator<Met>().deallocate(array_, NMAX);
}

monojet::MetCollection::reference
monojet::MetCollection::at(UInt_t _idx)
{
  if (_idx < size)
    return array_[_idx];
  throw std::out_of_range("MetCollection::at");
}

monojet::MetCollection::const_reference
monojet::MetCollection::at(UInt_t _idx) const
{
  if (_idx < size)
    return array_[_idx];
  throw std::out_of_range("MetCollection::at");
}

void
monojet::MetCollection::resize(UInt_t _size)
{
  if (size < NMAX) {
    size = _size;
    return;
  }
  throw std::length_error("MetCollection::resize");
}

void
monojet::MetCollection::setAddress(TTree& _tree)
{
  _tree.SetBranchAddress("met.size", &size);
  _tree.SetBranchAddress("met.met", met);
  _tree.SetBranchAddress("met.phi", phi);
}

void
monojet::MetCollection::book(TTree& _tree)
{
  _tree.Branch("met.size", &size, "size/i");
  _tree.Branch("met.met", met, "met[met.size]/F");
  _tree.Branch("met.phi", phi, "phi[met.size]/F");
}

monojet::ElectronCollection::ElectronCollection() :
  array_(std::allocator<Electron>().allocate(NMAX))
{
  for (unsigned iP(0); iP != NMAX; ++iP)
    new (array_ + iP) Electron(*this, iP);
}

monojet::ElectronCollection::~ElectronCollection()
{
  std::allocator<Electron>().deallocate(array_, NMAX);
}

monojet::ElectronCollection::reference
monojet::ElectronCollection::at(UInt_t _idx)
{
  if (_idx < size)
    return array_[_idx];
  throw std::out_of_range("ElectronCollection::at");
}

monojet::ElectronCollection::const_reference
monojet::ElectronCollection::at(UInt_t _idx) const
{
  if (_idx < size)
    return array_[_idx];
  throw std::out_of_range("ElectronCollection::at");
}

void
monojet::ElectronCollection::resize(UInt_t _size)
{
  if (size < NMAX) {
    size = _size;
    return;
  }
  throw std::length_error("ElectronCollection::resize");
}

void
monojet::ElectronCollection::setAddress(TTree& _tree)
{
  _tree.SetBranchAddress("electron.size", &size);
  _tree.SetBranchAddress("electron.pt", pt);
  _tree.SetBranchAddress("electron.eta", eta);
  _tree.SetBranchAddress("electron.phi", phi);
}

void
monojet::ElectronCollection::book(TTree& _tree)
{
  _tree.Branch("electron.size", &size, "size/i");
  _tree.Branch("electron.pt", pt, "pt[electron.size]/F");
  _tree.Branch("electron.eta", eta, "eta[electron.size]/F");
  _tree.Branch("electron.phi", phi, "phi[electron.size]/F");
}

monojet::MuonCollection::MuonCollection() :
  array_(std::allocator<Muon>().allocate(NMAX))
{
  for (unsigned iP(0); iP != NMAX; ++iP)
    new (array_ + iP) Muon(*this, iP);
}

monojet::MuonCollection::~MuonCollection()
{
  std::allocator<Muon>().deallocate(array_, NMAX);
}

monojet::MuonCollection::reference
monojet::MuonCollection::at(UInt_t _idx)
{
  if (_idx < size)
    return array_[_idx];
  throw std::out_of_range("MuonCollection::at");
}

monojet::MuonCollection::const_reference
monojet::MuonCollection::at(UInt_t _idx) const
{
  if (_idx < size)
    return array_[_idx];
  throw std::out_of_range("MuonCollection::at");
}

void
monojet::MuonCollection::resize(UInt_t _size)
{
  if (size < NMAX) {
    size = _size;
    return;
  }
  throw std::length_error("MuonCollection::resize");
}

void
monojet::MuonCollection::setAddress(TTree& _tree)
{
  _tree.SetBranchAddress("muon.size", &size);
  _tree.SetBranchAddress("muon.pt", pt);
  _tree.SetBranchAddress("muon.eta", eta);
  _tree.SetBranchAddress("muon.phi", phi);
}

void
monojet::MuonCollection::book(TTree& _tree)
{
  _tree.Branch("muon.size", &size, "size/i");
  _tree.Branch("muon.pt", pt, "pt[muon.size]/F");
  _tree.Branch("muon.eta", eta, "eta[muon.size]/F");
  _tree.Branch("muon.phi", phi, "phi[muon.size]/F");
}

monojet::TauCollection::TauCollection() :
  array_(std::allocator<Tau>().allocate(NMAX))
{
  for (unsigned iP(0); iP != NMAX; ++iP)
    new (array_ + iP) Tau(*this, iP);
}

monojet::TauCollection::~TauCollection()
{
  std::allocator<Tau>().deallocate(array_, NMAX);
}

monojet::TauCollection::reference
monojet::TauCollection::at(UInt_t _idx)
{
  if (_idx < size)
    return array_[_idx];
  throw std::out_of_range("TauCollection::at");
}

monojet::TauCollection::const_reference
monojet::TauCollection::at(UInt_t _idx) const
{
  if (_idx < size)
    return array_[_idx];
  throw std::out_of_range("TauCollection::at");
}

void
monojet::TauCollection::resize(UInt_t _size)
{
  if (size < NMAX) {
    size = _size;
    return;
  }
  throw std::length_error("TauCollection::resize");
}

void
monojet::TauCollection::setAddress(TTree& _tree)
{
  _tree.SetBranchAddress("tau.size", &size);
  _tree.SetBranchAddress("tau.pt", pt);
  _tree.SetBranchAddress("tau.eta", eta);
  _tree.SetBranchAddress("tau.phi", phi);
  _tree.SetBranchAddress("tau.disc", disc);
}

void
monojet::TauCollection::book(TTree& _tree)
{
  _tree.Branch("tau.size", &size, "size/i");
  _tree.Branch("tau.pt", pt, "pt[tau.size]/F");
  _tree.Branch("tau.eta", eta, "eta[tau.size]/F");
  _tree.Branch("tau.phi", phi, "phi[tau.size]/F");
  _tree.Branch("tau.disc", disc, "disc[tau.size]/F");
}

monojet::PhotonCollection::PhotonCollection() :
  array_(std::allocator<Photon>().allocate(NMAX))
{
  for (unsigned iP(0); iP != NMAX; ++iP)
    new (array_ + iP) Photon(*this, iP);
}

monojet::PhotonCollection::~PhotonCollection()
{
  std::allocator<Photon>().deallocate(array_, NMAX);
}

monojet::PhotonCollection::reference
monojet::PhotonCollection::at(UInt_t _idx)
{
  if (_idx < size)
    return array_[_idx];
  throw std::out_of_range("PhotonCollection::at");
}

monojet::PhotonCollection::const_reference
monojet::PhotonCollection::at(UInt_t _idx) const
{
  if (_idx < size)
    return array_[_idx];
  throw std::out_of_range("PhotonCollection::at");
}

void
monojet::PhotonCollection::resize(UInt_t _size)
{
  if (size < NMAX) {
    size = _size;
    return;
  }
  throw std::length_error("PhotonCollection::resize");
}

void
monojet::PhotonCollection::setAddress(TTree& _tree)
{
  _tree.SetBranchAddress("photon.size", &size);
  _tree.SetBranchAddress("photon.pt", pt);
  _tree.SetBranchAddress("photon.eta", eta);
  _tree.SetBranchAddress("photon.phi", phi);
}

void
monojet::PhotonCollection::book(TTree& _tree)
{
  _tree.Branch("photon.size", &size, "size/i");
  _tree.Branch("photon.pt", pt, "pt[photon.size]/F");
  _tree.Branch("photon.eta", eta, "eta[photon.size]/F");
  _tree.Branch("photon.phi", phi, "phi[photon.size]/F");
}


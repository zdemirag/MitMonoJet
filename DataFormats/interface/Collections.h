#ifndef monojet_Collections_h
#define monojet_Collections_h
#include "MitMonoJet/DataFormats/interface/Objects.h"

class TTree;

namespace monojet {

  class JetCollection {
  public:
    static UInt_t const NMAX = 32;
    typedef monojet::Jet value_type;
    typedef value_type& reference;
    typedef value_type const& const_reference;
    typedef monojet::Jet* iterator;
    typedef monojet::Jet const* const_iterator;

    JetCollection();
    ~JetCollection();

    reference at(UInt_t idx);
    const_reference at(UInt_t idx) const;
    void clear() { resize(0); }
    void resize(UInt_t size);
    iterator begin() { return array_; }
    const_iterator begin() const { return array_; }
    iterator end() { return array_ + size; }
    const_iterator end() const { return array_ + size; }

    void setAddress(TTree&);
    void book(TTree&);

    UInt_t size = 0;
    Float_t pt[NMAX];
    Float_t eta[NMAX];
    Float_t phi[NMAX];
    Float_t btag[NMAX];
    Float_t chFrac[NMAX];
    Float_t nhFrac[NMAX];
    Float_t neFrac[NMAX];

  private:
    value_type* array_;

  };

  class MetCollection {
  public:
    static UInt_t const NMAX = 1;
    typedef monojet::Met value_type;
    typedef value_type& reference;
    typedef value_type const& const_reference;
    typedef monojet::Met* iterator;
    typedef monojet::Met const* const_iterator;

    MetCollection();
    ~MetCollection();

    reference at(UInt_t idx);
    const_reference at(UInt_t idx) const;
    void clear() { resize(0); }
    void resize(UInt_t size);
    iterator begin() { return array_; }
    const_iterator begin() const { return array_; }
    iterator end() { return array_ + size; }
    const_iterator end() const { return array_ + size; }

    void setAddress(TTree&);
    void book(TTree&);

    UInt_t size = 0;
    Float_t met[NMAX];
    Float_t phi[NMAX];

  private:
    value_type* array_;

  };

  class ElectronCollection {
  public:
    static UInt_t const NMAX = 32;
    typedef monojet::Electron value_type;
    typedef value_type& reference;
    typedef value_type const& const_reference;
    typedef monojet::Electron* iterator;
    typedef monojet::Electron const* const_iterator;

    ElectronCollection();
    ~ElectronCollection();

    reference at(UInt_t idx);
    const_reference at(UInt_t idx) const;
    void clear() { resize(0); }
    void resize(UInt_t size);
    iterator begin() { return array_; }
    const_iterator begin() const { return array_; }
    iterator end() { return array_ + size; }
    const_iterator end() const { return array_ + size; }

    void setAddress(TTree&);
    void book(TTree&);

    UInt_t size = 0;
    Float_t pt[NMAX];
    Float_t eta[NMAX];
    Float_t phi[NMAX];

  private:
    value_type* array_;

  };

  class MuonCollection {
  public:
    static UInt_t const NMAX = 32;
    typedef monojet::Muon value_type;
    typedef value_type& reference;
    typedef value_type const& const_reference;
    typedef monojet::Muon* iterator;
    typedef monojet::Muon const* const_iterator;

    MuonCollection();
    ~MuonCollection();

    reference at(UInt_t idx);
    const_reference at(UInt_t idx) const;
    void clear() { resize(0); }
    void resize(UInt_t size);
    iterator begin() { return array_; }
    const_iterator begin() const { return array_; }
    iterator end() { return array_ + size; }
    const_iterator end() const { return array_ + size; }

    void setAddress(TTree&);
    void book(TTree&);

    UInt_t size = 0;
    Float_t pt[NMAX];
    Float_t eta[NMAX];
    Float_t phi[NMAX];

  private:
    value_type* array_;

  };

  class TauCollection {
  public:
    static UInt_t const NMAX = 32;
    typedef monojet::Tau value_type;
    typedef value_type& reference;
    typedef value_type const& const_reference;
    typedef monojet::Tau* iterator;
    typedef monojet::Tau const* const_iterator;

    TauCollection();
    ~TauCollection();

    reference at(UInt_t idx);
    const_reference at(UInt_t idx) const;
    void clear() { resize(0); }
    void resize(UInt_t size);
    iterator begin() { return array_; }
    const_iterator begin() const { return array_; }
    iterator end() { return array_ + size; }
    const_iterator end() const { return array_ + size; }

    void setAddress(TTree&);
    void book(TTree&);

    UInt_t size = 0;
    Float_t pt[NMAX];
    Float_t eta[NMAX];
    Float_t phi[NMAX];
    Float_t disc[NMAX];

  private:
    value_type* array_;

  };

  class PhotonCollection {
  public:
    static UInt_t const NMAX = 32;
    typedef monojet::Photon value_type;
    typedef value_type& reference;
    typedef value_type const& const_reference;
    typedef monojet::Photon* iterator;
    typedef monojet::Photon const* const_iterator;

    PhotonCollection();
    ~PhotonCollection();

    reference at(UInt_t idx);
    const_reference at(UInt_t idx) const;
    void clear() { resize(0); }
    void resize(UInt_t size);
    iterator begin() { return array_; }
    const_iterator begin() const { return array_; }
    iterator end() { return array_ + size; }
    const_iterator end() const { return array_ + size; }

    void setAddress(TTree&);
    void book(TTree&);

    UInt_t size = 0;
    Float_t pt[NMAX];
    Float_t eta[NMAX];
    Float_t phi[NMAX];

  private:
    value_type* array_;

  };

}

#endif

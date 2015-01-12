// $Id: XlFatJet.cc,v 1.1 2008/06/04 09:08:36 loizides Exp $

#include "MitMonoJet/DataTree/interface/XlFatJet.h"

ClassImp(mithep::XlFatJet)

namespace mithep {

const XlSubJet* XlFatJet::SubJet(UInt_t i, XlSubJet::ESubJetType t) const {
    UInt_t NJets = fSubJets.Entries();
    UInt_t counter = 0;
    for (UInt_t j=0;j<NJets;j++) {
        if (fSubJets.At(j)->SubJetType() == t) {
            counter++;
            if (counter == i) {
                return fSubJets.At(j);
            }
        }
    }
    return NULL;
}

UInt_t XlFatJet::NSubJets(XlSubJet::ESubJetType t) const {
    UInt_t NJets = fSubJets.Entries();
    UInt_t counter = 0;
    for (UInt_t j=0;j<NJets;j++) {
        if (fSubJets.At(j)->SubJetType() == t) {
            counter++;
        }
    }
    return counter;
}

UInt_t XlFatJet::NTopSubJets() const {  return NSubJets(XlSubJet::ESubJetType::eTop); }

UInt_t XlFatJet::NVSubJets() const {  return NSubJets(XlSubJet::ESubJetType::eV); }

}
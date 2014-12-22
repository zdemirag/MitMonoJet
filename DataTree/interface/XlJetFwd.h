// $Id:$

#ifndef MITMONOJET_DATATREE_XLJETFWD_H
#define MITMONOJET_DATATREE_XLJETFWD_H

#include "MitAna/DataCont/interface/CollectionFwd.h"
#include "MitAna/DataCont/interface/ArrayFwd.h"
#include "MitAna/DataCont/interface/ObjArrayFwd.h"

namespace mithep {
  class XlJet;
  typedef Collection<XlJet> XlJetCol;
  typedef Array<XlJet> XlJetArr;
  typedef ObjArray<XlJet> XlJetOArr;
}
#endif

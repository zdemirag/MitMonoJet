// $Id:$

#ifndef MITMONOJET_DATATREE_XLMETFWD_H
#define MITMONOJET_DATATREE_XLMETFWD_H

#include "MitAna/DataCont/interface/CollectionFwd.h"
#include "MitAna/DataCont/interface/ArrayFwd.h"
#include "MitAna/DataCont/interface/ObjArrayFwd.h"

namespace mithep {
  class XlMet;
  typedef Collection<XlMet> XlMetCol;
  typedef Array<XlMet> XlMetArr;
  typedef ObjArray<XlMet> XlMetOArr;
}
#endif

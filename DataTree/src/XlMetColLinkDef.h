// $Id: XlMetColLinkDef.h,v 1.1 2009/06/15 15:00:15 loizides Exp $

#ifndef MITMONOJET_DATATREE_XLMETCOLLINKDEF_H
#define MITMONOJET_DATATREE_XLMETCOLLINKDEF_H

#include "MitAna/DataCont/interface/Ref.h"
#include "MitMonoJet/DataTree/interface/XlMetCol.h"
#endif

#ifndef __CINT__
# define _R__UNIQUEIDENTIFIER_ XlMetCol
# define _R__JOIN3_(F,X,Y) _NAME3_(F,X,Y)
# undef _R__UNIQUE_
# define _R__UNIQUE_(X) _R__JOIN3_( _R__UNIQUEIDENTIFIER_,X,__LINE__)
#endif

#ifdef __CINT__
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;
#pragma link C++ nestedclass;
#pragma link C++ nestedtypedef;
#pragma link C++ namespace mithep;

#pragma link C++ class mithep::XlMet+;
#pragma link C++ class mithep::Collection<mithep::XlMet>+;
#pragma link C++ class mithep::Array<mithep::XlMet>+;
#pragma link C++ class mithep::ObjArray<mithep::XlMet>+;
#pragma link C++ class mithep::Ref<mithep::XlMet>+;
#pragma link C++ typedef mithep::XlMetCol;
#pragma link C++ typedef mithep::XlMetArr;
#pragma link C++ typedef mithep::XlMetOArr;
#endif

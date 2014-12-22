// $Id: XlJetColLinkDef.h,v 1.1 2009/06/15 15:00:15 loizides Exp $

#ifndef MITMONOJET_DATATREE_XLJETCOLLINKDEF_H
#define MITMONOJET_DATATREE_XLJETCOLLINKDEF_H

#include "MitAna/DataCont/interface/Ref.h"
#include "MitMonoJet/DataTree/interface/XlJetCol.h"
#endif

#ifndef __CINT__
# define _R__UNIQUEIDENTIFIER_ XlJetCol
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

#pragma link C++ class mithep::XlJet+;
#pragma link C++ class mithep::Collection<mithep::XlJet>+;
#pragma link C++ class mithep::Array<mithep::XlJet>+;
#pragma link C++ class mithep::ObjArray<mithep::XlJet>+;
#pragma link C++ class mithep::Ref<mithep::XlJet>+;
#pragma link C++ typedef mithep::XlJetCol;
#pragma link C++ typedef mithep::XlJetArr;
#pragma link C++ typedef mithep::XlJetOArr;
#endif

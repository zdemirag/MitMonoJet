// $Id: XsIsoParticleColLinkDef.h,v 1.1 2009/06/15 15:00:15 loizides Exp $

#ifndef MITMONOJET_DATATREE_XSISOPARTICLECOLLINKDEF_H
#define MITMONOJET_DATATREE_XSISOPARTICLECOLLINKDEF_H

#include "MitAna/DataCont/interface/RefArray.h"
#include "MitAna/DataCont/interface/Ref.h"
#include "MitMonoJet/DataTree/interface/XsIsoParticleCol.h"
#endif

#ifndef __CINT__
# define _R__UNIQUEIDENTIFIER_ XsIsoParticleCol
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

#pragma link C++ class mithep::XsIsoParticle+;
#pragma link C++ class mithep::Collection<mithep::XsIsoParticle>+;
#pragma link C++ class mithep::Array<mithep::XsIsoParticle>+;
#pragma link C++ class mithep::ObjArray<mithep::XsIsoParticle>+;
#pragma link C++ class mithep::Ref<mithep::XsIsoParticle>+;
#pragma link C++ class mithep::RefArray<mithep::XsIsoParticle>+;
#pragma link C++ typedef mithep::XsIsoParticleCol;
#pragma link C++ typedef mithep::XsIsoParticleArr;
#pragma link C++ typedef mithep::XsIsoParticleOArr;
#endif

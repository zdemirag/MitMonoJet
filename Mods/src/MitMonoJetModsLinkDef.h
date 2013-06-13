// $Id: MitMonoPhotonModsLinkDef.h,v 1.1 2013/05/31 18:20:56 dimatteo Exp $

#ifndef MITMONOJET_MODS_LINKDEF_H
#define MITMONOJET_MODS_LINKDEF_H
#include "MitMonoJet/Mods/interface/MonoJetTreeWriter.h"
#endif

#ifdef __CINT__
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;
#pragma link C++ nestedclass;
#pragma link C++ nestedtypedef;
#pragma link C++ namespace mithep;

#pragma link C++ class mithep::MonoJetTreeWriter+;
#pragma link C++ class mithep::MonoJetEvent+;
#endif

// $Id: MitMonoJetModsLinkDef.h,v 1.2 2013/07/11 14:28:16 dimatteo Exp $

#ifndef MITMONOJET_MODS_LINKDEF_H
#define MITMONOJET_MODS_LINKDEF_H
#include "MitMonoJet/Mods/interface/MonoJetTreeWriter.h"
#include "MitMonoJet/Mods/interface/HLTEvtSelMod.h"
#endif

#ifdef __CINT__
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;
#pragma link C++ nestedclass;
#pragma link C++ nestedtypedef;
#pragma link C++ namespace mithep;

#pragma link C++ class mithep::MonoJetTreeWriter+;
#pragma link C++ class mithep::HLTEvtSelMod+;
#pragma link C++ class mithep::MonoJetEventHLT+;
#endif

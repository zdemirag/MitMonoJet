// $Id: MitMonoJetModsLinkDef.h,v 1.3 2013/07/19 05:54:39 dimatteo Exp $

#ifndef MITMONOJET_MODS_LINKDEF_H
#define MITMONOJET_MODS_LINKDEF_H
#include "MitMonoJet/Mods/interface/MonoJetTreeWriter.h"
#include "MitMonoJet/Mods/interface/HltEvtSelMod.h"
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
#pragma link C++ class mithep::HltEvtSelMod+;
#pragma link C++ class mithep::MonoJetEventHlt+;

#pragma link C++ class mithep::HLTEvtSelMod+;
#pragma link C++ class mithep::MonoJetEventHLT+;

#endif

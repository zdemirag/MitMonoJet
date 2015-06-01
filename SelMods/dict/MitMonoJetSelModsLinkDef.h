// $Id: MitMonoJetSelModsLinkDef.h,v 1.1 2013/06/13 03:31:18 dimatteo Exp $

#ifndef MITMONOJET_SELMODS_LINKDEF_H
#define MITMONOJET_SELMODS_LINKDEF_H

#include "MitMonoJet/SelMods/interface/BoostedVAnalysisMod.h"
#include "MitMonoJet/SelMods/interface/MonoJetAnalysisMod.h"
#endif
 
#ifdef __CINT__
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;
#pragma link C++ nestedclass;
#pragma link C++ nestedtypedef;
#pragma link C++ namespace mithep;

#pragma link C++ class mithep::BoostedVAnalysisMod+;
#pragma link C++ class mithep::ObjArray<mithep::TriggerMask>+;
#pragma link C++ class mithep::MonoJetAnalysisMod+;
#endif

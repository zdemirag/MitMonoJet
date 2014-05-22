
#ifndef MITMONOJET_UTILS_LINKDEF_H
#define MITMONOJET_UTILS_LINKDEF_H

#include "MitMonoJet/Utils/interface/Qjets.h"
#include "MitMonoJet/Utils/interface/QjetsPlugin.h"
#endif

#ifdef __CINT__
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;
#pragma link C++ nestedclass;
#pragma link C++ nestedtypedef;

#pragma link C++ class JetDistanceCompare;
#pragma link C++ class Qjets;
#pragma link C++ class QjetsPlugin;
#endif

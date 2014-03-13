//--------------------------------------------------------------------------------------------------
// Perform a plot task using a specified set of samples. Nice and clean.
//
// Authors: C.Paus                                                                        (Aug 2010)
//--------------------------------------------------------------------------------------------------
#include "MitPlots/Style/interface/MitStyle.h"
#include "MitPlots/Plot/interface/Plot.h"

using namespace std;
using namespace mithep;

//==================================================================================================
void simplePlots(double lumi = 4500.0)
{
  // setup graphics stuff before starting
  MitStyle::Init();

  // plot from TTree
  TString nTuple   = "MJetTree";
  TString variable = "genJet1Tau2/genJet1Tau1";
  TString cut      = TString("genJet1Pt>300 && genJet1M>65. && genJet1M<95. &&") +
                     TString("genJet1Eta<1.3 && genJet1Eta>-1.3");

  printf("\n Looking at tree (%s)\n   plotting variable (%s)\n   with cuts: (%s)\n\n",
	 nTuple.Data(),variable.Data(),cut.Data());

  plot(nTuple.Data(),"genTau2OverTau1","#tau_{2}/#tau_{1}",0,0.07,1.0,0.,0.,1,lumi,variable,cut,93);

  return;
}

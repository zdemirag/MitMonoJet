#include "MitMonoJet/DataFormats/interface/Event.h"
#include "TTree.h"

void
monojet::Event::setAddress(TTree& _tree)
{
  _tree.SetBranchAddress("runNumber", &runNumber);
  _tree.SetBranchAddress("lumiNumber", &lumiNumber);
  _tree.SetBranchAddress("eventNumber", &eventNumber);

  jets.setAddress(_tree);
  mets.setAddress(_tree);
  electrons.setAddress(_tree);
  muons.setAddress(_tree);
  taus.setAddress(_tree);
  photons.setAddress(_tree);
}

void
monojet::Event::book(TTree& _tree)
{
  _tree.Branch("runNumber", &runNumber, "runNumber/i");
  _tree.Branch("lumiNumber", &lumiNumber, "lumiNumber/i");
  _tree.Branch("eventNumber", &eventNumber, "eventNumber/i");

  jets.book(_tree);
  mets.book(_tree);
  electrons.book(_tree);
  muons.book(_tree);
  taus.book(_tree);
  photons.book(_tree);
}


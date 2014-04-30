//--------------------------------------------------------------------------------------------------
// $Id: FillerXlJets.h,v 1.9 2011/03/01 17:27:22 mzanetti Exp $
//
// FillerXlJets
//
// This module process a collection of input jet, compute the substrucure
// and fill two output collections of XlFatJets and relative XlSubJets
//
// Authors: L.DiMatteo
//--------------------------------------------------------------------------------------------------

#ifndef MITMONOJET_TREEFILLER_FILLERXLJETS_H
#define MITMONOJET_TREEFILLER_FILLERXLJETS_H

#include "fastjet/PseudoJet.hh"
#include "fastjet/JetDefinition.hh"
#include "fastjet/GhostedAreaSpec.hh"
#include "fastjet/AreaDefinition.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/tools/Pruner.hh"

#include "fastjet/contrib/Njettiness.hh"
#include "fastjet/contrib/Nsubjettiness.hh"
#include "fastjet/contrib/NjettinessPlugin.hh"
#include "MitMonoJet/DataTree/interface/XlFatJetFwd.h"
#include "MitMonoJet/DataTree/interface/XlFatJet.h"
#include "MitMonoJet/DataTree/interface/XlSubJetFwd.h"

#include "MitAna/TreeMod/interface/BaseMod.h"
#include "MitAna/DataTree/interface/JetCol.h"
#include "MitAna/DataTree/interface/PFCandidateCol.h"

namespace mithep
{
  class FillerXlJets : public BaseMod
  {
    public:
      FillerXlJets(const char *name = "FillerXlJets",
                   const char *title = "XlJets Filler module");
      ~FillerXlJets();

      void SetJetsName(const char *n)      { fJetsName = n;       }
      void SetJetsFromBranch(bool b)       { fJetsFromBranch = b; }
 
      void FillXlFatJets(std::vector<fastjet::PseudoJet> &fjFatJets);
      void FillXlSubJets(std::vector<fastjet::PseudoJet> &fjSubJets, XlFatJet *pFatJet);

    protected:
      void Process();
      void SlaveBegin();
      void SlaveTerminate();
  
    private:
      Bool_t fIsData;                      //is this data or MC?
      Bool_t fBTaggingActive;              //=true if BTagging info is filled
      Bool_t fQGTaggingActive;             //=true if QGTagging info is filled
      Bool_t fPublishOutput;               //=true if output collection are published

      TString fJetsName;                   //(i) name of input jets
      Bool_t fJetsFromBranch;              //are input jets from Branch?
      const JetCol *fJets;                 //input jets

      TString fPfCandidatesName;           //(i) name of PF candidates coll
      Bool_t fPfCandidatesFromBranch;      //are PF candidates from Branch?
      const PFCandidateCol *fPfCandidates; //particle flow candidates coll handle
 
      TString XlFatJetsName;               //name of output XlFatJets collection
      XlFatJetArr *XlFatJets;              //array of XlFatJets
      TString XlSubJetsName;               //name of output XlSubJets collection
      XlSubJetArr *XlSubJets;              //array of XlSubJets
      
      // Objects from fastjet we want to use
      double fConeSize;
      int fPrune;                      //apply pruning: 0-no, 1-standard CMS
      fastjet::Pruner *fPruner;
      fastjet::JetDefinition *fCAJetDef;
      fastjet::GhostedAreaSpec *fActiveArea;
      fastjet::AreaDefinition *fAreaDefinition;

      ClassDef(FillerXlJets, 0)        //XlJets, Fat and Sub, filler      
  };
}
#endif

//--------------------------------------------------------------------------------------------------
// $Id: FillerXlFatJets.h,v 1.9 2011/03/01 17:27:22 mzanetti Exp $
//
// FillerXlFatJets
//
// This module process a collection of input jet, compute the substrucure
// and fill two output collections of fXlFatJets and relative fXlSubJets
//
// Authors: L.DiMatteo
//--------------------------------------------------------------------------------------------------

#ifndef MITMONOJET_TREEFILLER_FILLERXLFATJETS_H
#define MITMONOJET_TREEFILLER_FILLERXLFATJETS_H

#include <TVector2.h>

#include "fastjet/PseudoJet.hh"
#include "fastjet/JetDefinition.hh"
#include "fastjet/GhostedAreaSpec.hh"
#include "fastjet/AreaDefinition.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/tools/Pruner.hh"
#include "fastjet/tools/Filter.hh"

#include "fastjet/contrib/Njettiness.hh"
#include "fastjet/contrib/EnergyCorrelator.hh"
#include "fastjet/contrib/Nsubjettiness.hh"
#include "fastjet/contrib/NjettinessPlugin.hh"
#include "fastjet/contrib/SoftDrop.hh"
#include "MitMonoJet/DataTree/interface/XlFatJetFwd.h"
#include "MitMonoJet/DataTree/interface/XlFatJet.h"
#include "MitMonoJet/DataTree/interface/XlSubJetFwd.h"
#include "MitMonoJet/DataTree/interface/XlSubJet.h"

#include "MitAna/TreeMod/interface/BaseMod.h"
#include "MitAna/DataTree/interface/JetCol.h"
#include "MitAna/DataTree/interface/PFCandidateCol.h"
#include "MitAna/DataTree/interface/PFJet.h"
#include "MitAna/DataTree/interface/PileupEnergyDensityCol.h"
#include "MitAna/DataTree/interface/VertexCol.h"
#include "MitPhysics/Utils/interface/QGTagger.h"

namespace mithep
{
  class FillerXlFatJets : public BaseMod
  {
    public:
      FillerXlFatJets(const char *name = "FillerXlFatJets",
                   const char *title = "XlFatJets Filler module");
      ~FillerXlFatJets();

      void IsData(Bool_t b)                { fIsData = b;           }
      void FillVSubJets(Bool_t b)          { fFillVSubJets = b;     }
      void FillTopSubJets(Bool_t b)        { fFillTopSubJets = b;   }
      void NSubDeclustering(Bool_t b)      { fNSubDeclustering = b; }
      void SetBtaggingOn(Bool_t b)         { fBTaggingActive = b;   }
      void SetfQGTaggingOn(Bool_t b)       { fQGTaggingActive = b;  }
      void setTopTaggingOn(Bool_t b)       { fTopTaggingActive = b; }
      void SetQGTaggerCHS(bool b)          { fQGTaggerCHS = b;      }
      void PublishOutput(Bool_t b)         { fPublishOutput = b;    }
                                                                    
      void SetProcessNJets(UInt_t n)       { fProcessNJets = n;     } 
                                                                    
      void SetJetsName(const char *n)      { fJetsName = n;         }
      void SetJetsFromBranch(Bool_t b)     { fJetsFromBranch = b;   }
                                                                    
      void SetFatJetsName(const char *n)   { fXlFatJetsName = n;    }
      void SetSubJetsName(const char *n)   { fXlSubJetsName = n;    }
                                                                    
      void SetSoftDropZCut(double d)       { fSoftDropZCut = d;     }
      void SetSoftDropR0(double d)      { fSoftDropR0 = d;    }
      void SetPruneZCut(double d)          { fPruneZCut = d;        }
      void SetPruneDistCut(double d)       { fPruneDistCut = d;     }
      void SetFilterN(int n)               { fFilterN = n;          }
      void SetFilterRad(double d)          { fFilterRad = d;        }
      void SetTrimRad(double d)            { fTrimRad = d;          }
      void SetTrimPtFrac(double d)         { fTrimPtFrac = d;       }
      void SetConeSize(double d)           { fConeSize = d;         }

 
    protected:
      void Process();
      void SlaveBegin();
      void SlaveTerminate();
 
      void FillXlFatJet (const PFJet *pPFJet);
      void FillXlSubJets(std::vector<fastjet::PseudoJet> &fjSubJets,
                         XlFatJet *pFatJet,XlSubJet::ESubJetType t);
      // Jet collection helpers
      std::vector <fastjet::PseudoJet>   
            Sorted_by_pt_min_pt(std::vector <fastjet::PseudoJet> &jets,  
                                float jetPtMin);      
      // QJets volatility helpers
      void   GetJetConstituents(fastjet::PseudoJet &jet, std::vector <fastjet::PseudoJet> &constits,  
                                float constitsPtMin);
      double GetQjetVolatility (std::vector<fastjet::PseudoJet> &constits, int QJetsN = 25, int seed = 12345);
      double FindRMS           (std::vector<float> qjetmasses);
      double FindMean          (std::vector<float> qjetmasses); 
      // Subjet QGTagging helpers
      void   FillSubjetQGTagging(fastjet::PseudoJet &jet, float constitsPtMin, 
                                 XlSubJet *pSubJet, XlFatJet *pFatJet);
      // Color pull helpers
      TVector2 GetPull(fastjet::PseudoJet &jet, float constitsPtMin);
      double GetPullAngle(std::vector<fastjet::PseudoJet> &fjSubJets, float constitsPtMin);

    private:
      Bool_t fIsData;                      //is this data or MC?
      Bool_t fFillVSubJets;                //=true if V-subjets are stored (2-prom structure)
      Bool_t fFillTopSubJets;              //=true if top-subjets are stored (3-prom structure)
      Bool_t fNSubDeclustering;            //=true if subjets declustering via n-subjettiness
                                           //=false if subjets declustering via pruned CA (CMS standard)
      Bool_t fBTaggingActive;              //=true if BTagging info is filled
      Bool_t fQGTaggingActive;             //=true if QGTagging info is filled
      Bool_t fTopTaggingActive;            //=true if CMSTopTagger info is filled
      Bool_t fQGTaggerCHS;                 //=true if QGTagging weights are taken from CHS
      Bool_t fPublishOutput;               //=true if output collection are published

      UInt_t  fProcessNJets;               //number of input jets processed by fastjet

      TString fJetsName;                   //(i) name of input jets
      Bool_t fJetsFromBranch;              //are input jets from Branch?
      const JetCol *fJets;                 //input jets

      TString fPfCandidatesName;           //(i) name of PF candidates coll
      Bool_t fPfCandidatesFromBranch;      //are PF candidates from Branch?
      const PFCandidateCol *fPfCandidates; //particle flow candidates coll handle

      TString fPileUpDenName;              //(i) name of PU energy density coll
      Bool_t fPileUpDenFromBranch;         //is PU energy density from Branch?
      const PileupEnergyDensityCol *fPileUpDen; //PU energy density coll handle

      TString fVertexesName;               //(i) name of vertex coll
      Bool_t fVertexesFromBranch;          //are vertexex from Branch?
      const VertexCol *fVertexes;          //vertex coll handle
 
      TString fXlFatJetsName;              //name of output fXlFatJets collection
      XlFatJetArr *fXlFatJets;             //array of fXlFatJets
      TString fXlSubJetsName;              //name of output fXlSubJets collection
      XlSubJetArr *fXlSubJets;             //array of fXlSubJets
      
      // Objects from fastjet we want to use
      fastjet::Pruner *fPruner;
      fastjet::Filter *fFilterer;
      fastjet::Filter *fTrimmer;           //no this is not a typo trimmers belong to fastjet Filter class
      double fSoftDropZCut;                //soft-drop Z cut
      double fSoftDropR0;               //soft-drop angular distance normalisation
      double fPruneZCut;                   //pruning Z cut
      double fPruneDistCut;                //pruning distance cut 
      int fFilterN;                        //number of subjets after filtering
      double fFilterRad;                   //filtered subjets radius
      double fTrimRad;                     //trimmed subjet radius
      double fTrimPtFrac;                  //trimmed subjet pt fraction
      double fConeSize;                    //fastjet clustering radius
      fastjet::JetDefinition *fCAJetDef;   //fastjet clustering definition
      fastjet::GhostedAreaSpec *fActiveArea;
      fastjet::AreaDefinition *fAreaDefinition;
      
      // QG tagger 
      QGTagger *fQGTagger;                 //QGTagger calculator
      
      // Counters : used to initialize seed for QJets volatility
      Long64_t fCounter;

      ClassDef(FillerXlFatJets, 0)         //XlJets, Fat and Sub, filler      
  };
}
#endif

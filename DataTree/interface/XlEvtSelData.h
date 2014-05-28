//--------------------------------------------------------------------------------------------------
// $Id: EvtSelData.h,v $
//
// EvtSelData
//
// Class to store event properties on which one may want to select clean events. 
// This extension also store information about the preselection module decision.
//
// Authors: C.Loizides, L.Di Matteo
//--------------------------------------------------------------------------------------------------

#ifndef MITMONOJET_DATATREE_XLEVTSELDATA_H
#define MITMONOJET_DATATREE_XLEVTSELDATA_H

#include <TString.h> 
#include "MitAna/DataTree/interface/EvtSelData.h"
#include "MitAna/DataTree/interface/DataObject.h"

namespace mithep 
{
  class XlEvtSelData : public EvtSelData
  {
    public:
      XlEvtSelData(const char *name=0) : 
                fPreselWord(0), fName(name) {}                 
      XlEvtSelData(const EvtSelData & e, const char *name=0) : 
                EvtSelData(e),
                fPreselWord(0), fName(name) {}                 

      Int_t                 preselWord()               const { return fPreselWord; }
      const char*           GetName()                  const;
 
      void                  SetPreselWord(Int_t t)           { fPreselWord  = t;   }
      void                  SetName(const char* name)        { fName = name;       }

    protected:
      Int_t                 fPreselWord;   //word containing the bit-values of the preselection filters
      TString               fName;         //name of the object (should be specified in DataBase)

    ClassDef(XlEvtSelData, 1) // XlEvtSelData class
  };
}

//--------------------------------------------------------------------------------------------------
inline const char* mithep::XlEvtSelData::GetName() const
{
   // Return name of this collection.
   // if no name, return the collection class name.
   if (fName.Length() > 0) return fName.Data();
   return ClassName();
}

#endif



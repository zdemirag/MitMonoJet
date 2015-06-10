#include <TFile.h>
#include <TLorentzVector.h>
#include "TMVA/Tools.h"

#include "MitAna/DataTree/interface/MCParticleCol.h"
#include "MitMonoJet/Utils/interface/DiJetMVA.h"

using namespace mithep;

//--------------------------------------------------------------------------------------------------
DiJetMVA::DiJetMVA() :
  fLowPtMethodName ("BDTG"),
  fHighPtMethodName("BDTG"),
  fIsInitialized(kFALSE),
  fApplyQGCorrection(kTRUE),
  fDiJetPtThr(160.),
  fLowPtReader(0),
  fHighPtReader(0),
  fQGSyst(0)
{      
}

//--------------------------------------------------------------------------------------------------
DiJetMVA::~DiJetMVA() 
{
  delete fLowPtReader;
  delete fHighPtReader;
  delete fQGSyst;
}

//--------------------------------------------------------------------------------------------------
void DiJetMVA::Initialize( TString lowPtWeights,
			                     TString highPtWeights, 
			                     std::string QGCorrDatabase ) 
{
  fIsInitialized = kTRUE;
  
  // If requested setup QGCorrector
  if (fApplyQGCorrection) {
    fQGSyst = new QGSyst;
    fQGSyst->ReadDatabaseDoubleMin(QGCorrDatabase);
  }  
  // Set up low pt MVA 
  fLowPtReader = new TMVA::Reader( "!Color:!Silent:Error" );  
  fLowPtReader->AddSpectator("isSig",   &fMVAVar_isSig);
  fLowPtReader->AddVariable ("ptmjj",   &fMVAVar_ptmjj);
  fLowPtReader->AddVariable ("pull1ang",&fMVAVar_pull1ang);
  fLowPtReader->AddVariable ("pull2ang",&fMVAVar_pull2ang);
  fLowPtReader->AddVariable ("qgid1",   &fMVAVar_qgid1);
  fLowPtReader->AddVariable ("qgid2",   &fMVAVar_qgid2);
  fLowPtReader->AddVariable ("mdrop",   &fMVAVar_mdrop);
  fLowPtReader->BookMVA(fLowPtMethodName, lowPtWeights );

  // Set up high pt MVA 
  fHighPtReader = new TMVA::Reader( "!Color:!Silent:Error" );  
  fHighPtReader->AddSpectator("isSig",   &fMVAVar_isSig);
  fHighPtReader->AddVariable ("ptmjj",   &fMVAVar_ptmjj);
  fHighPtReader->AddVariable ("pull1ang",&fMVAVar_pull1ang);
  fHighPtReader->AddVariable ("pull2ang",&fMVAVar_pull2ang);
  fHighPtReader->AddVariable ("qgid1",   &fMVAVar_qgid1);
  fHighPtReader->AddVariable ("qgid2",   &fMVAVar_qgid2);
  fHighPtReader->AddVariable ("mdrop",   &fMVAVar_mdrop);
  fHighPtReader->BookMVA(fHighPtMethodName, highPtWeights );

  // Say what we are doing   
  std::cout << "DiJet V-Tag MVA Initialization\n";
  std::cout << "MethodName : " << fLowPtMethodName << std::endl;
  std::cout << "Weight file: " << lowPtWeights << std::endl;
  
  return;
}

//--------------------------------------------------------------------------------------------------
Double_t DiJetMVA::MVAValue( const XlJet *jet1, const XlJet *jet2,
                             const Float_t rho, 
                             const MCParticleCol *genParticles, 
                             Bool_t printDebug, float* vars )
//--------------------------------------------------------------------------------------------------
{  
  if(!fIsInitialized) { 
    std::cout << "Error: DiJetMVA not properly initialized.\n"; 
    return -1.;
  }

  // Prepare the vectors
  TLorentzVector vJet1,vJet2;
  vJet1.SetPtEtaPhiE(jet1->Pt(),jet1->Eta(),jet1->Phi(),jet1->E());
  vJet2.SetPtEtaPhiE(jet2->Pt(),jet2->Eta(),jet2->Phi(),jet2->E());
  TLorentzVector vDijet = vJet1 + vJet2;
   
  // Prepare the variables with dummy spectator
  fMVAVar_ptmjj    = vDijet.Pt()/vDijet.M();
  fMVAVar_pull1ang = computePullAngle(jet1,jet2);
  fMVAVar_pull2ang = computePullAngle(jet2,jet1);
  fMVAVar_qgid1    = jet1->QGTag();
  fMVAVar_qgid2    = jet2->QGTag();
  fMVAVar_mdrop    = TMath::Max(vJet1.M(), vJet2.M())/vDijet.M()*(vJet1.DeltaR(vJet2));
  // Adjust QG if needed and MC is valid
  if (fApplyQGCorrection && genParticles) {
    int iflav1 = JetPartonMatch(jet1,genParticles,0.5);  
    std::string sflav1="";
    if(abs(iflav1) >= 1 && abs(iflav1) <= 5) { sflav1 = "quark"; } 
    if(abs(iflav1) == 21)                    { sflav1 = "gluon"; }
    if (abs(iflav1) >= 1 && (abs(iflav1) <= 5 || abs(iflav1) == 21) )
      fMVAVar_qgid1 = fQGSyst->Smear(vJet1.Pt(), vJet1.Eta(), rho, fMVAVar_qgid1, sflav1);

    int iflav2 = JetPartonMatch(jet2,genParticles,0.5);  
    std::string sflav2="";
    if(abs(iflav2) >= 1 && abs(iflav2) <= 5) { sflav2 = "quark"; } 
    if(abs(iflav2) == 21)                    { sflav2 = "gluon"; }
    if (abs(iflav2) >= 1 && (abs(iflav2) <= 5 || abs(iflav2) == 21) )
      fMVAVar_qgid2 = fQGSyst->Smear(vJet2.Pt(), vJet2.Eta(), rho, fMVAVar_qgid2, sflav2);
  }

  // Evaluate MVA
  Double_t mvaval = (vDijet.Pt()<fDiJetPtThr) ? fHighPtReader->EvaluateMVA(fHighPtMethodName) : fHighPtReader->EvaluateMVA(fHighPtMethodName);

  // Debug if needed
  if (printDebug) {
    std::cout << "Debug DiJet MVA: "
              << jet1->Pt() << " " << jet1->Eta() << " " << jet1->Phi() << " , "
              << jet2->Pt() << " " << jet2->Eta() << " " << jet2->Phi() << " : "
              << fMVAVar_ptmjj            << " " 
              << fMVAVar_pull1ang         << " " 
              << fMVAVar_pull2ang         << " " 
              << fMVAVar_qgid1            << " " 
              << fMVAVar_qgid2            << " "  
              << fMVAVar_mdrop            << " " 
              << " === : === "
              << mvaval 
              << std::endl;
  }

  // Store variables if needed
  if (vars) {
    vars[0] = fMVAVar_ptmjj   ;
    vars[1] = fMVAVar_pull1ang;
    vars[2] = fMVAVar_pull2ang;
    vars[3] = fMVAVar_qgid1   ;
    vars[4] = fMVAVar_qgid2   ; 
    vars[5] = fMVAVar_mdrop   ;
  }
  
  return mvaval;
}

//--------------------------------------------------------------------------------------------------
Float_t DiJetMVA::computePullAngle(const XlJet *basisJet, const XlJet *testJet)
{

  // change coordinates of jet axes for jet pull angle computation
  TVector2 jet2_in_jet1coord(testJet->Rapidity() - basisJet->Rapidity(), MathUtils::DeltaPhi(testJet,basisJet));

  TVector2 pullvec;
  pullvec.Set(basisJet->PullY(), basisJet->PullPhi());

  return pullvec.DeltaPhi(jet2_in_jet1coord);
}

//--------------------------------------------------------------------------------------------------
Int_t DiJetMVA::JetPartonMatch(const XlJet *jet,
                               const MCParticleCol *genParticles,
                               Float_t deltaR)
{
  float minDr = 999.;
  unsigned int pId = 0;
      
  // Loop on all stable MC particles and perform the matching
  for (UInt_t i=0; i<genParticles->GetEntries(); ++i) {
    const MCParticle *p = genParticles->At(i);
    if (p->Status()!=3)
      continue;
    // compute the dR    
    float thisDr = MathUtils::DeltaR(*jet, *p);
    // standard check for any parton    
    if (thisDr < minDr) {
      minDr = thisDr;
      pId = p->AbsPdgId();
    }
  }

  // return matched parton only if minDr less that user dr limit
  if (minDr < deltaR)
    return pId;
  else 
    return 0;  
}

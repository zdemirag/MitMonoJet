#include <TROOT.h>
#include "TMath.h"
#include "Math/VectorUtil.h"
#include "Math/LorentzVector.h"
typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > LorentzVector; 

void functions(){
}

float getMass(float e1, float e2, float phi1, float phi2){
  return TMath::Sqrt(2.*e1*e2*(1-cos(phi1-phi2)));
}

float getMass(float e1, float e2, float phi1, float phi2, float p1, float p2){
  return TMath::Sqrt(e1*e1-p1*p1 + e2*e2-p2*p2 + 2.*(e1*e2-p1*p2*cos(phi1-phi2)));
}

float deltaPhi(float phi1, float phi2) {
  float PHI = TMath::Abs(phi1-phi2);
  if (PHI<=TMath::Pi())
    return PHI;
  else
    return 2*TMath::Pi()-PHI;
}

float deltaR(float phi1, float eta1, float phi2, float eta2) {
  return TMath::Sqrt((eta2-eta1)*(eta2-eta1)+deltaPhi(phi1,phi2)*deltaPhi(phi1,phi2));
}

float vectorSumPhi(float px1, float py1, float px2, float py2){
  float phi = atan((py1+py2)/(px1+px2));
  if ((px1+px2)>0) return phi;
  else if ((py1+py2)>0) return phi + TMath::Pi();
  else return phi - TMath::Pi();
}

float vectorSumPt(float pt1, float phi1, float pt2, float phi2){
  return sqrt( pow(pt1*cos(phi1) + pt2*cos(phi2),2) +
	       pow(pt1*sin(phi1) + pt2*sin(phi2),2) );
}

float vectorSum3Pt(float pt1, float phi1, float pt2, float phi2,float pt3, float phi3){
  return sqrt( pow(pt1*cos(phi1) + pt2*cos(phi2) + pt3*cos(phi3),2) +
	       pow(pt1*sin(phi1) + pt2*sin(phi2) + pt3*sin(phi3),2) );
}

float vectorSumMass(float px1, float py1, float pz1, float px2, float py2, float pz2) {
  double E1 = TMath::Sqrt(px1*px1 + py1*py1 + pz1*pz1);
  double E2 = TMath::Sqrt(px2*px2 + py2*py2 + pz2*pz2);
  double cosTheta = (px1*px2 + py1*py2 + pz1*pz2)/ (E1*E2);
  return TMath::Sqrt(2*E1*E2*(1-cosTheta));
}

float transverseMass(float lepPt, float lepPhi, float met,  float metPhi) {
  double cosDPhi = cos(deltaPhi(lepPhi,metPhi));
  return sqrt(2*lepPt*met*(1-cosDPhi));
}

float caloMet1l(float pt, float phi, float met, float metPhi){
  return sqrt( pow(pt*cos(phi) + met*cos(metPhi),2) +
	       pow(pt*sin(phi) + met*sin(metPhi),2));
}

float caloMet2l(float pt1, float phi1, float pt2, float phi2, float met, float metPhi){
  return sqrt( pow(pt1*cos(phi1) + pt2*cos(phi2) + met*cos(metPhi),2) +
	       pow(pt1*sin(phi1) + pt2*sin(phi2) + met*sin(metPhi),2));
}

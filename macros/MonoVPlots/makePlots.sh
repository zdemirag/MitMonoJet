#!/bin/bash
#---------------------------------------------------------------------------------------------------
# Prepare weighted histogram for nice plotting
#
#                                                                         L.Di Matteo (Nov 14, 2014)
#---------------------------------------------------------------------------------------------------
if [ $1 = "sig" ];
then
  # signal region plots
  root -b -q ../rootlogon_monojet.C makeRawPlot.C+'(19700,1,"metRaw",4,55,250,800)'
  root -b -q ../rootlogon_monojet.C makeRawPlot.C+'(19700,1,"fjet1.M()",4,40,0,200)'
  root -b -q ../rootlogon_monojet.C makeRawPlot.C+'(19700,1,"fjet1.Pt()",4,55,250,800)'
  root -b -q ../rootlogon_monojet.C makeRawPlot.C+'(19700,1,"bdt_all",5,40,-1.,1.)'
elif [ $1 = "zll" ];
then
  # zll cr region plots, HORRIBLE fix for weird name
  root -b -q ../rootlogon_monojet.C makeRawPlot.C+'(19700,2,"met",2,55,250,800)'
  mv BDT_Zll_TMath\:\:Sqrt\(TMath\:\:Power\(metRaw\*TMath\:\:Cos\(metRawPhi\ +\ lep1.Px\(\ +\ lep2.Px\(\,2\ +\ TMath\:\:Power\(metRaw\*TMath\:\:Sin\(metRawPhi\ +\ lep1.Py\(\ +\ lep2.Py\(\,2.root BDT_Zll_met.root
  root -b -q ../rootlogon_monojet.C makeRawPlot.C+'(19700,2,"fjet1.M()",2,40,0,200)'
  root -b -q ../rootlogon_monojet.C makeRawPlot.C+'(19700,2,"fjet1.Pt()",2,55,250,800)'
  root -b -q ../rootlogon_monojet.C makeRawPlot.C+'(19700,2,"bdt_all",6,40,-1.,1.)'
  exit 0
fi
  
exit 0

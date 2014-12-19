#!/bin/bash
#---------------------------------------------------------------------------------------------------
# Prepare weighted histogram for nice plotting
#
#                                                                         L.Di Matteo (Nov 14, 2014)
#---------------------------------------------------------------------------------------------------
# For each region following the same schema
# 1) boostedV plots
# 2) relevant vtag variables
# 3) baseline plots
# 4) inclusive plots (counting only)
if [ $1 = "sig" ]  || [ $1 = "all" ];
then
  # signal region plots
  root -b -q ../rootlogon_monojet.C makeRawPlot.C+'(19700,1,"metRaw",5,55,250,1000,1)'
  root -b -q ../rootlogon_monojet.C makeRawPlot.C+'(19700,1,"fjet1.Pt()",5,55,250,1000,1)'
  root -b -q ../rootlogon_monojet.C makeRawPlot.C+'(19700,1,"isData",5,2,0,2)'

  root -b -q ../rootlogon_monojet.C makeRawPlot.C+'(19700,1,"fjet1MassPruned",14,40,0,200)'

  root -b -q ../rootlogon_monojet.C makeRawPlot.C+'(19700,1,"metRaw",10,55,200,1000,1)'
  root -b -q ../rootlogon_monojet.C makeRawPlot.C+'(19700,1,"jet1.Pt()",10,55,150,1000,1)'
  root -b -q ../rootlogon_monojet.C makeRawPlot.C+'(19700,1,"fjet1Tau2/fjet1Tau1",10,40,0.,1.)'
  root -b -q ../rootlogon_monojet.C makeRawPlot.C+'(19700,1,"fjet1MassPruned",10,40,0,200)'
  root -b -q ../rootlogon_monojet.C makeRawPlot.C+'(19700,1,"isData",10,2,0,2)'

  root -b -q ../rootlogon_monojet.C makeRawPlot.C+'(19700,1,"isData",6,2,0,2)'
fi
if [ $1 = "zll" ] || [ $1 = "all" ];
then
  # zll cr region plots, HORRIBLE fix for weird name
  root -b -q ../rootlogon_monojet.C makeRawPlot.C+'(19700,2,"met",2,35,250,1000,1)'
  root -b -q ../rootlogon_monojet.C makeRawPlot.C+'(19700,2,"fjet1.Pt()",2,35,250,1000,1)'
  root -b -q ../rootlogon_monojet.C makeRawPlot.C+'(19700,2,"isData",2,2,0,2)'
  
  root -b -q ../rootlogon_monojet.C makeRawPlot.C+'(19700,2,"fjet1MassPruned",15,40,0,200)'

  root -b -q ../rootlogon_monojet.C makeRawPlot.C+'(19700,2,"met",11,55,200,1000,1)'
  root -b -q ../rootlogon_monojet.C makeRawPlot.C+'(19700,2,"jet1.Pt()",11,55,150,1000,1)'
  root -b -q ../rootlogon_monojet.C makeRawPlot.C+'(19700,2,"fjet1Tau2/fjet1Tau1",11,40,0.,1.)'
  root -b -q ../rootlogon_monojet.C makeRawPlot.C+'(19700,2,"fjet1MassPruned",11,40,0,200)'
  root -b -q ../rootlogon_monojet.C makeRawPlot.C+'(19700,2,"isData",11,2,0,2)'

  root -b -q ../rootlogon_monojet.C makeRawPlot.C+'(19700,2,"isData",7,2,0,2)'
fi
if [ $1 = "wlv" ] || [ $1 = "all" ];
then
  # wlv cr region plots, HORRIBLE fix for weird name
  root -b -q ../rootlogon_monojet.C makeRawPlot.C+'(19700,3,"met",3,35,250,1000,1)'
  root -b -q ../rootlogon_monojet.C makeRawPlot.C+'(19700,3,"fjet1.Pt()",3,45,250,1000,1)'
  root -b -q ../rootlogon_monojet.C makeRawPlot.C+'(19700,3,"isData",3,2,0,2)'

  root -b -q ../rootlogon_monojet.C makeRawPlot.C+'(19700,3,"fjet1MassPruned",16,40,0,200)'

  root -b -q ../rootlogon_monojet.C makeRawPlot.C+'(19700,3,"met",12,55,200,1000,1)'
  root -b -q ../rootlogon_monojet.C makeRawPlot.C+'(19700,3,"jet1.Pt()",12,55,150,1000,1)'
  root -b -q ../rootlogon_monojet.C makeRawPlot.C+'(19700,3,"fjet1Tau2/fjet1Tau1",12,40,0.,1.)'
  root -b -q ../rootlogon_monojet.C makeRawPlot.C+'(19700,3,"fjet1MassPruned",12,40,0,200)'
  root -b -q ../rootlogon_monojet.C makeRawPlot.C+'(19700,3,"isData",12,2,0,2)'

  root -b -q ../rootlogon_monojet.C makeRawPlot.C+'(19700,3,"isData",8,2,0,2)'
fi
if [ $1 = "pj" ] || [ $1 = "all" ];
then
  # pj cr region plots, HORRIBLE fix for weird name
  root -b -q ../rootlogon_monojet.C makeRawPlot.C+'(19700,4,"met",4,55,250,1000,1)'
  root -b -q ../rootlogon_monojet.C makeRawPlot.C+'(19700,4,"fjet1.Pt()",4,55,250,1000,1)'
  root -b -q ../rootlogon_monojet.C makeRawPlot.C+'(19700,4,"isData",4,40,0,200)'

  root -b -q ../rootlogon_monojet.C makeRawPlot.C+'(19700,4,"fjet1MassPruned",17,40,0,200)'
  
  root -b -q ../rootlogon_monojet.C makeRawPlot.C+'(19700,4,"met",13,75,200,1000,1)'
  root -b -q ../rootlogon_monojet.C makeRawPlot.C+'(19700,4,"jet1.Pt()",13,75,200,1000,1)'
  root -b -q ../rootlogon_monojet.C makeRawPlot.C+'(19700,4,"fjet1Tau2/fjet1Tau1",13,40,0.,1.)'
  root -b -q ../rootlogon_monojet.C makeRawPlot.C+'(19700,4,"fjet1MassPruned",13,40,0,200)'
  root -b -q ../rootlogon_monojet.C makeRawPlot.C+'(19700,4,"isData",13,40,0,200)'

  root -b -q ../rootlogon_monojet.C makeRawPlot.C+'(19700,4,"isData",9,40,0,200)'
fi
  
exit 0

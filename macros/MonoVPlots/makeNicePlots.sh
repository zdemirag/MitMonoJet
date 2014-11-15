#!/bin/bash
#---------------------------------------------------------------------------------------------------
# Make the plots with CMS style with all bells and whistles 
#
#                                                                         L.Di Matteo (Nov 14, 2014)
#---------------------------------------------------------------------------------------------------
if [ $1 = "sig" ];
then
  # signal
  root -l -q -b finalPlotMonoJet.C+'(0,1,"MET ","GeV","BDT_Met_metRaw.root","met_sig",1,19.7,1,0)'
  root -l -q -b finalPlotMonoJet.C+'(0,1,"Jet Mass ","GeV","BDT_Met_fjet1.M().root","mass_sig",0,19.7,1,0)'
  root -l -q -b finalPlotMonoJet.C+'(0,1,"Jet P_{T} ","GeV","BDT_Met_fjet1.Pt().root","pt_sig",1,19.7,1,0)'
  root -l -q -b finalPlotMonoJet.C+'(0,1,"V-Tag MVA ","","BDT_Met_bdt_all.root","bdt_sig",1,19.7,1,0)'
elif [ $1 = "zll" ];
then
  # zll CR region plotss
  root -l -q -b finalPlotMonoJet.C+'(1,1,"MET ","GeV","BDT_Zll_met.root","met_zll",1,19.7,0,0)'
  root -l -q -b finalPlotMonoJet.C+'(1,1,"Jet Mass ","GeV","BDT_Zll_fjet1.M().root","mass_zll",0,19.7,0,0)'
  root -l -q -b finalPlotMonoJet.C+'(1,1,"Jet P_{T} ","GeV","BDT_Zll_fjet1.Pt().root","pt_zll",1,19.7,0,0)'
  root -l -q -b finalPlotMonoJet.C+'(1,1,"V-Tag MVA ","","BDT_Zll_bdt_all.root","bdt_zll",1,19.7,0,0)'
  exit 0
fi
  
exit 0
  

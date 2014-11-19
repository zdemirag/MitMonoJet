#!/bin/bash
#---------------------------------------------------------------------------------------------------
# Make the plots with CMS style with all bells and whistles 
#
#                                                                         L.Di Matteo (Nov 14, 2014)
#---------------------------------------------------------------------------------------------------
if [ $1 = "sig" ] || [ $1 = "all" ];
then
  # signal
  root -l -q -b finalPlotMonoJet.C+'(0,1,"MET ","GeV","BDT_Met_metRaw.root","met_sig",1,19.7,1,0)'
  root -l -q -b finalPlotMonoJet.C+'(0,1,"Jet Mass ","GeV","BDT_Met_fjet1.M().root","mass_sig",0,19.7,1,0)'
  root -l -q -b finalPlotMonoJet.C+'(0,1,"Jet P_{T} ","GeV","BDT_Met_fjet1.Pt().root","pt_sig",1,19.7,1,0)'
  root -l -q -b finalPlotMonoJet.C+'(0,1,"V-Tag MVA ","","BDT_Met_noVTagCut_bdt_all.root","bdt_sig_novtag",1,19.7,1,0)'
  root -l -q -b finalPlotMonoJet.C+'(0,1,"MET ","GeV","BDT_Met_noVTagCut_metRaw.root","met_sig_novtag",1,19.7,1,0)'
  root -l -q -b finalPlotMonoJet.C+'(0,1,"Jet Mass ","GeV","BDT_Met_noVTagCut_fjet1.M().root","mass_sig_novtag",0,19.7,1,0)'
  root -l -q -b finalPlotMonoJet.C+'(0,1,"Jet P_{T} ","GeV","BDT_Met_noVTagCut_fjet1.Pt().root","pt_novtag",1,19.7,1,0)'
fi
if [ $1 = "zll" ] || [ $1 = "all" ];
then
  # zll CR region plots
  root -l -q -b finalPlotMonoJet.C+'(1,1,"MET ","GeV","BDT_Zll_met.root","met_zll",1,19.7,0,0)'
  root -l -q -b finalPlotMonoJet.C+'(1,1,"Jet Mass ","GeV","BDT_Zll_fjet1.M().root","mass_zll",0,19.7,0,0)'
  root -l -q -b finalPlotMonoJet.C+'(1,1,"Jet P_{T} ","GeV","BDT_Zll_fjet1.Pt().root","pt_zll",1,19.7,0,0)'
  root -l -q -b finalPlotMonoJet.C+'(1,1,"V-Tag MVA ","","BDT_Zll_noVTagCut_bdt_all.root","bdt_zll_novtag",1,19.7,0,0)'
  root -l -q -b finalPlotMonoJet.C+'(1,1,"MET ","GeV","BDT_Zll_noVTagCut_met.root","met_zll_novtag",1,19.7,0,0)'
  root -l -q -b finalPlotMonoJet.C+'(1,1,"Jet Mass ","GeV","BDT_Zll_noVTagCut_fjet1.M().root","mass_zll_novtag",0,19.7,0,0)'
  root -l -q -b finalPlotMonoJet.C+'(1,1,"Jet P_{T} ","GeV","BDT_Zll_noVTagCut_fjet1.Pt().root","pt_zll_novtag",1,19.7,0,0)'
fi
if [ $1 = "wlv" ] || [ $1 = "all" ];
then
  # wlv CR region plots
  root -l -q -b finalPlotMonoJet.C+'(2,1,"MET ","GeV","BDT_Wlv_met.root","met_wlv",1,19.7,0,0)'
  root -l -q -b finalPlotMonoJet.C+'(2,1,"Jet Mass ","GeV","BDT_Wlv_fjet1.M().root","mass_wlv",0,19.7,0,0)'
  root -l -q -b finalPlotMonoJet.C+'(2,1,"Jet P_{T} ","GeV","BDT_Wlv_fjet1.Pt().root","pt_wlv",1,19.7,0,0)'
  root -l -q -b finalPlotMonoJet.C+'(2,1,"V-Tag MVA ","","BDT_Wlv_noVTagCut_bdt_all.root","bdt_wlv_novtag",1,19.7,0,0)'
  root -l -q -b finalPlotMonoJet.C+'(2,1,"MET ","GeV","BDT_Wlv_noVTagCut_met.root","met_wlv_novtag",1,19.7,0,0)'
  root -l -q -b finalPlotMonoJet.C+'(2,1,"Jet Mass ","GeV","BDT_noVTagCut_Wlv_fjet1.M().root","mass_wlv_novtag",0,19.7,0,0)'
  root -l -q -b finalPlotMonoJet.C+'(2,1,"Jet P_{T} ","GeV","BDT_Wlv_noVTagCut_fjet1.Pt().root","pt_wlv_novtag",1,19.7,0,0)'
fi
if [ $1 = "pj" ] || [ $1 = "all" ];
then
  # pj CR region plots
  root -l -q -b finalPlotMonoJet.C+'(3,1,"MET ","GeV","BDT_Pj_met.root","met_pj",1,19.7,0,1)'
  root -l -q -b finalPlotMonoJet.C+'(3,1,"Jet Mass ","GeV","BDT_Pj_fjet1.M().root","mass_pj",0,19.7,0,1)'
  root -l -q -b finalPlotMonoJet.C+'(3,1,"Jet P_{T} ","GeV","BDT_Pj_fjet1.Pt().root","pt_pj",1,19.7,0,1)'
  root -l -q -b finalPlotMonoJet.C+'(3,1,"V-Tag MVA ","","BDT_Pj_noVTagCut_bdt_all.root","BDT_Pj_novtag",1,19.7,0,1)'
  root -l -q -b finalPlotMonoJet.C+'(3,1,"MET ","GeV","BDT_Pj_noVTagCut_met.root","met_pj_novtag",1,19.7,0,1)'
  root -l -q -b finalPlotMonoJet.C+'(3,1,"Jet Mass ","GeV","BDT_noVTagCut_pj_fjet1.M().root","mass_pj_novtag",0,19.7,0,1)'
  root -l -q -b finalPlotMonoJet.C+'(3,1,"Jet P_{T} ","GeV","BDT_Pj_noVTagCut_fjet1.Pt().root","pt_pj_novtag",1,19.7,0,1)'
fi
  
exit 0
  

#!/bin/bash
#---------------------------------------------------------------------------------------------------
# Make the plots with CMS style with all bells and whistles 
#
#                                                                         L.Di Matteo (Nov 14, 2014)
#---------------------------------------------------------------------------------------------------
if [ $1 = "sig" ] || [ $1 = "all" ];
then
  # signal
  root -l -q -b finalPlotMonoJet.C+'(0,1,"MET ","GeV","BDT_Met_metRaw.root","met_sig",1,19.7,1,0,1)'
  root -l -q -b finalPlotMonoJet.C+'(0,1,"C-A 0.8 Jet P_{T} ","GeV","BDT_Met_fjet1.Pt().root","pt_sig",1,19.7,1,0,1)'
  
  root -l -q -b finalPlotMonoJet.C+'(0,1,"Jet Mass_{prun} ","GeV","BDT_Met_nomasscut_fjet1MassPruned.root","mass_sig_nomasscut",0,19.7,1,0)'
  
  root -l -q -b finalPlotMonoJet.C+'(0,1,"MET ","GeV","BDT_Met_baseline_metRaw.root","met_sig_baseline",1,19.7,1,0,1)'
  root -l -q -b finalPlotMonoJet.C+'(0,1,"Ak_{T} 0.5 Jet P_{T} ","GeV","BDT_Met_baseline_jet1.Pt().root","pt_sig_baseline",1,19.7,1,0,1)'
  root -l -q -b finalPlotMonoJet.C+'(0,1,"C-A 0.8 Jet #tau_{2}/#tau_{1} ","","BDT_Met_baseline_fjet1Tau2DIVfjet1Tau1.root","t2t1_sig_baseline",0,19.7,1,0)'
  root -l -q -b finalPlotMonoJet.C+'(0,1,"Jet Mass_{prun} ","GeV","BDT_Met_baseline_fjet1MassPruned.root","mass_sig_baseline",0,19.7,1,0)'  
fi
if [ $1 = "zll" ] || [ $1 = "all" ];
then
  # zll CR region plots
  root -l -q -b finalPlotMonoJet.C+'(1,1,"MET ","GeV","BDT_Zll_met.root","met_zll",1,19.7,0,1,1)'
  root -l -q -b finalPlotMonoJet.C+'(1,1,"C-A 0.8 Jet P_{T} ","GeV","BDT_Zll_fjet1.Pt().root","pt_zll",1,19.7,0,1,1)'

  root -l -q -b finalPlotMonoJet.C+'(1,1,"Jet Mass_{prun} ","GeV","BDT_Zll_nomasscut_fjet1MassPruned.root","mass_zll_nomasscut",0,19.7,0,1)'

  root -l -q -b finalPlotMonoJet.C+'(1,1,"MET ","GeV","BDT_Zll_baseline_met.root","met_zll_baseline",1,19.7,0,1,1)'
  root -l -q -b finalPlotMonoJet.C+'(1,1,"Ak_{T} 0.5 Jet P_{T} ","GeV","BDT_Zll_baseline_jet1.Pt().root","pt_zll_baseline",1,19.7,0,1,1)'
  root -l -q -b finalPlotMonoJet.C+'(1,1,"C-A 0.8 Jet #tau_{2}/#tau_{1} ","","BDT_Zll_baseline_fjet1Tau2DIVfjet1Tau1.root","t2t1_zll_baseline",0,19.7,0,1)'
  root -l -q -b finalPlotMonoJet.C+'(1,1,"Jet Mass_{prun} ","GeV","BDT_Zll_baseline_fjet1MassPruned.root","mass_zll_baseline",0,19.7,0,1)'
fi
if [ $1 = "wlv" ] || [ $1 = "all" ];
then
  # wlv CR region plots
  root -l -q -b finalPlotMonoJet.C+'(2,1,"MET ","GeV","BDT_Wlv_met.root","met_wlv",1,19.7,0,1,1)'
  root -l -q -b finalPlotMonoJet.C+'(2,1,"C-A 0.8 Jet P_{T} ","GeV","BDT_Wlv_fjet1.Pt().root","pt_wlv",1,19.7,0,1,1)'

  root -l -q -b finalPlotMonoJet.C+'(2,1,"Jet Mass_{prun} ","GeV","BDT_Wlv_nomasscut_fjet1MassPruned.root","mass_wlv_nomasscut",0,19.7,0,1)'

  root -l -q -b finalPlotMonoJet.C+'(2,1,"MET ","GeV","BDT_Wlv_baseline_met.root","met_wlv_baseline",1,19.7,0,1,1)'
  root -l -q -b finalPlotMonoJet.C+'(2,1,"Ak_{T} 0.5 Jet P_{T} ","GeV","BDT_Wlv_baseline_jet1.Pt().root","pt_wlv_baseline",1,19.7,0,1,1)'
  root -l -q -b finalPlotMonoJet.C+'(2,1,"C-A 0.8 Jet #tau_{2}/#tau_{1} ","","BDT_Wlv_baseline_fjet1Tau2DIVfjet1Tau1.root","t2t1_wlv_baseline",0,19.7,0,1)'
  root -l -q -b finalPlotMonoJet.C+'(2,1,"Jet Mass_{prun} ","GeV","BDT_Wlv_baseline_fjet1MassPruned.root","mass_wlv_baseline",0,19.7,0,1)'
fi
if [ $1 = "pj" ] || [ $1 = "all" ];
then
  # pj CR region plots
  root -l -q -b finalPlotMonoJet.C+'(3,1,"MET ","GeV","BDT_Pj_met.root","met_pj",1,19.7,0,1,1)'
  root -l -q -b finalPlotMonoJet.C+'(3,1,"C-A 0.8 Jet P_{T} ","GeV","BDT_Pj_fjet1.Pt().root","pt_pj",1,19.7,0,1,1)'

  root -l -q -b finalPlotMonoJet.C+'(3,1,"Jet Mass_{prun} ","GeV","BDT_Pj_nomasscut_fjet1MassPruned.root","mass_pj_nomasscut",0,19.7,0,1)'

  root -l -q -b finalPlotMonoJet.C+'(3,1,"MET ","GeV","BDT_Pj_baseline_met.root","met_pj_baseline",1,19.7,0,1,1)'
  root -l -q -b finalPlotMonoJet.C+'(3,1,"Ak_{T} 0.5 Jet P_{T} ","GeV","BDT_Pj_baseline_jet1.Pt().root","pt_pj_baseline",1,19.7,0,1,1)'  
  root -l -q -b finalPlotMonoJet.C+'(3,1,"C-A 0.8 Jet #tau_{2}/#tau_{1} ","","BDT_Pj_baseline_fjet1Tau2DIVfjet1Tau1.root","t2t1_pj_baseline",0,19.7,0,1)'
  root -l -q -b finalPlotMonoJet.C+'(3,1,"Jet Mass_{prun} ","GeV","BDT_Pj_baseline_fjet1MassPruned.root","mass_pj_baseline",0,19.7,0,1)'
fi
  
exit 0
  

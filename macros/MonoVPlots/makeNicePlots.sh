# signal
root -l -q -b finalPlotMonoJet.C+'(0,1,"MET ","GeV","BDT_Met_metRaw.root","met_sig",1,19.7,1,0)'
root -l -q -b finalPlotMonoJet.C+'(0,1,"Jet Mass ","GeV","BDT_Met_fjet1.M().root","mass_sig",0,19.7,1,0)'
root -l -q -b finalPlotMonoJet.C+'(0,1,"Jet P_{T} ","GeV","BDT_Met_fjet1.Pt().root","pt_sig",1,19.7,1,0)'
root -l -q -b finalPlotMonoJet.C+'(0,1,"V-Tag MVA ","","BDT_Met_bdt_all.root","bdt_sig",1,19.7,1,0)'

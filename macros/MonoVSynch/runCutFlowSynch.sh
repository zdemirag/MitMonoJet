#!/bin/bash
#---------------------------------------------------------------------------------------------------
# Produce synchronization files for cutflow comparison
#
#                                                                            L.Di Matteo (Dec, 2014)
#---------------------------------------------------------------------------------------------------
# TTbar
# - prepare the Bavanti ntuples
root -b -q ../rootlogon_monojet.C ../runBavantiBoostedV_Synch.C+'("0000","noskim","s12-h125inv-gf-v7c","t2mit/filefi/032","/home/cmsprod/catalog","boostedv",3000)'
root -b -q ../rootlogon_monojet.C ../runBavantiBoostedV_Synch.C+'("0000","noskim","GJets_HT-400ToInf_8TeV-madgraph+Summer12_DR53X-PU_S10_START53_V7A-v1+AODSIM","t2mit/filefi/032","/home/cmsprod/catalog","boostedv",300)'
# - prepare the input file for the flat step
ls *s12-h125inv-gf-v7c*noskim*0000*ntuple* > inputBavanti.txt
root -b -q ../rootlogon_monojet.C ../runFlatBoostedV.C+
mv boostedv-v9_s12-ttj-v1-v7a_noskim_flatntuple.root ./s12-h125inv-gf-v7c.root

ls *GJets_HT-400ToInf_8TeV-mad*noskim*0000*ntuple* > inputBavanti.txt
root -b -q ../rootlogon_monojet.C ../runFlatBoostedV.C+
mv boostedv-v9_s12-ttj-v1-v7a_noskim_flatntuple.root ./s12-gjets.root
# - create the synch tables
./cutflowDMSTree.py -f s12-h125inv-gf-v7c.root
./cutflowDMSTree.py -f s12-gjets.root -s gjets

# - check the min and max run
#./listeventsDMSTree.py -f s12-h125inv-gf-v7a.root -n -1 > runs.txt
# Cleanup
#rm ./*.root

exit 0

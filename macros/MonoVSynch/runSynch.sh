#!/bin/bash
#---------------------------------------------------------------------------------------------------
# Produce synchronization files
#
#                                                                            L.Di Matteo (Dec, 2014)
#---------------------------------------------------------------------------------------------------
# TTbar
# - prepare the Bavanti ntuples
root -b -q ../rootlogon_monojet.C runBavantiBoostedV_Synch.C+'("0000","noskim","s12-ttj-sl-v1-v7c","t2mit/filefi/032","/home/cmsprod/catalog","boostedv",10)'
# - prepare the input file for the flat step
ls *ttj-sl*ntuple* > inputBavanti.txt
# - run the flat step
root -b -q ../rootlogon_monojet.C ../runFlatBoostedV.C+
# - create the synch file
./dumpDMSTree.py -f DMSTreeWriter_tmp.root > BambuDump-tt.txt

# Gamma+jets
# - prepare the Bavanti ntuples
root -b -q ../rootlogon_monojet.C runBavantiBoostedV_Synch.C+'("0000","noskim","GJets_HT-400ToInf_8TeV-madgraph+Summer12_DR53X-PU_S10_START53_V7A-v1+AODSIM","t2mit/filefi/032","/home/cmsprod/catalog","boostedv",10)'
# - prepare the input file for the flat step
ls *GJets_HT-400*ntuple* > inputBavanti.txt
# - run the flat step
root -b -q ../rootlogon_monojet.C ../runFlatBoostedV.C+
# - create the synch file
./dumpDMSTree.py -f DMSTreeWriter_tmp.root > BambuDump-gjets.txt

# Cleanup
rm ./*.root

exit 0

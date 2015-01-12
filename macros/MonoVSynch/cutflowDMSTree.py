#! /usr/bin/env python

#----------------------------------------------------------------------
# Sample usage
# ./cutflowDMSTree.py -f /scratch4/dimatteo/cms/hist/boostedv-v10/merged-test/boostedv-v10_WH_WToQQ_HToInv_M-125_8TeV_pythia6+Summer12_DR53X-PU_S10_START53_V19-v1+AODSIM_noskim_flatntuple.root
#----------------------------------------------------------------------

import os, re, array, ROOT
import sys
from optparse import OptionParser

# - S U B R O U T I N E S -----------------------------------------------------------------------------

# Make table header
def makeTableHeader(tot_entries):
  print ' {0:<40s} & {1:>20s} \\\ \hline\hline'\
  .format('cut description','selected events') 

  print ' {0:<40s} & {1:20d} \\\ \hline'\
  .format('all events',tot_entries) 

  return

# Make table raw
def makeTableRaw(tree,histo,cut_string,cut_name):
  tree.Draw('isData >> ' + histo.GetName(),cut_string,'goff')
  n_sel = int(histo.GetSumOfWeights())
  print ' {0:<40s} & {1:20d} \\\ \hline'\
  .format(cut_name,n_sel) 

# - M A I N ----------------------------------------------------------------------------------------

# Prepare the command line parser
parser = OptionParser()
parser.add_option("-f", "--file", dest="input_file", default='boostedv_s12-sch-v7a_noskim_flatntuple.root',
                  help="input root file [default: %default]")
parser.add_option("-t", "--treename", dest="input_tree", default='DMSTree',
                  help="root tree name [default: %default]")
(options, args) = parser.parse_args()

# Get all the root classes
from ROOT import *

# Open the correct input file and get the event tree
input_file = TFile.Open(options.input_file)
if input_file:
  print 'INFO - Opening input root file: ' + options.input_file
  
else:
  print 'ERROR - Cannot open input root file: ' + options.input_file + ' , exiting!'
  raise SystemExit
  
input_tree = input_file.FindObjectAny(options.input_tree)
if input_tree:
  print 'INFO - Opening root tree: ' + options.input_tree
  
else:
  print 'ERROR - Cannot open root tree: ' + options.input_tree + ' , exiting!'
  raise SystemExit

# Check the total number of entries (BAMBU only!)
n_tot_entries = int(input_file.FindObjectAny("hDAllEvents").GetEntries())
print 'INFO - Total entries: ' + str(n_tot_entries)

# Check the number of entries in the tree
n_entries = input_tree.GetEntriesFast()
print 'INFO - Input tree entries: ' + str(n_entries)

# Prepare cutflow table header
makeTableHeader(n_tot_entries)

# Prepare histogram for cutflow printout
h_counts = TH1F( 'h_counts', 'Counting histogram based on isData', 2, 0, 2)
# Prepare cut string
cut_string = ''

# Preselection cut (BAMBU only!)
cut_string = cut_string + '(preselWord & (1<<3))'
makeTableRaw(input_tree,h_counts,cut_string,'HLT and rawPFMET > 150 and BAMBUpresel')

# Jet pt cut
cut_string = cut_string + '&&(jet1.Pt() > 150)'
makeTableRaw(input_tree,h_counts,cut_string,'AK5 pt > 150')

# Jet eta cut
cut_string = cut_string + '&&(abs(jet1.eta()) < 2.5)'
makeTableRaw(input_tree,h_counts,cut_string,'abs(AK5_{eta}) < 2.5')

# MET cut
cut_string = cut_string + '&&(metRaw > 200)'
makeTableRaw(input_tree,h_counts,cut_string,'rawPFMET > 200')

# Njets cut
cut_string = cut_string + '&&(njets < 3)'
makeTableRaw(input_tree,h_counts,cut_string,'nAK5 < 3')

# Delta phi jet1,jet2
cut_string = cut_string + '&&(njets == 1 || TMath::ACos(TMath::Cos(jet2.phi() - jet1.phi())) < 2.0)'
makeTableRaw(input_tree,h_counts,cut_string,'dPhi(AK5_1,AK5_2) < 2.0')

# Lepton veto
cut_string = cut_string + '&&(nlep == 0)'
makeTableRaw(input_tree,h_counts,cut_string,'nlep == 0')

# Photon veto
cut_string = cut_string + '&&(nphotons == 0)'
makeTableRaw(input_tree,h_counts,cut_string,'nphotons == 0')

# Tau veto
cut_string = cut_string + '&&(ntaus == 0)'
makeTableRaw(input_tree,h_counts,cut_string,'ntaus == 0')

# Met noise filters
cut_string = cut_string + '&&(metFiltersWord == 511 || metFiltersWord == 1023)'
makeTableRaw(input_tree,h_counts,cut_string,'MET noise filters')

# Jet noise filters
cut_string = cut_string + '&&(jet1CHF > 0.2 && jet1NHF < 0.7 && jet1NEMF < 0.7)'
makeTableRaw(input_tree,h_counts,cut_string,'AK5_1 noise filters')


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
  print ' {0:<40s} & {1:>20s} & {2:>20s} & {3:>20s} \\\ \hline\hline'\
  .format('cut description','selected events','absolute efficiency','relative efficiency') 

  print ' {0:<40s} & {1:20d} & {2:20.2f} & {3:20.2f} \\\ \hline'\
  .format('all events',tot_entries,100,100) 

  return

# Make table raw
def makeTableRaw(tree,histo,cut_string,cut_name,tot_entries,n_sel_previous):
  tree.Draw('isData >> ' + histo.GetName(),cut_string,'goff')
  n_sel = int(histo.GetSumOfWeights())
  print ' {0:<40s} & {1:20d} & {2:20.2f} & {3:20.2f} \\\ \hline'\
  .format(cut_name,n_sel,float(n_sel)/float(tot_entries)*100.,float(n_sel)/float(n_sel_previous)*100.) 
  
  histo.Reset()
  
  return n_sel

# - M A I N ----------------------------------------------------------------------------------------

# Prepare the command line parser
parser = OptionParser()
parser.add_option("-f", "--file", dest="input_file", default='boostedv_s12-sch-v7a_noskim_flatntuple.root',
                  help="input root file [default: %default]")
parser.add_option("-t", "--treename", dest="input_tree", default='DMSTree',
                  help="root tree name [default: %default]")
parser.add_option("-s", "--selection", dest="input_selection", default='signal',
                  help="selection type [default: %default]")
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
print '<verbatim>'
makeTableHeader(n_tot_entries)

# Prepare histogram for cutflow printout
h_counts = TH1F( 'h_counts', 'Counting histogram based on isData', 2, 0, 2)
# Prepare cut string
cut_string = ''
# Prepare counter
counter = n_tot_entries

# HLT cut (BAMBU only!)
if (options.input_selection == 'signal'):
  cut_string = cut_string + '((trigger & (1<<0)) > 0 || (trigger & (1<<1)) > 0)'
else:
  cut_string = cut_string + '((trigger & (1<<3)) > 0)'
  
counter = makeTableRaw(input_tree,h_counts,cut_string,'HLT',n_tot_entries,counter)

# Jet pt cut
cut_string = cut_string + '&&(jet1.Pt() > 150)'
counter = makeTableRaw(input_tree,h_counts,cut_string,'AK5 pt > 150',n_tot_entries,counter)

# MET cut
if (options.input_selection == 'signal'):
  cut_string = cut_string + '&&(metRaw > 200)'
else:
  cut_string = cut_string + '&&(met > 200)'
counter = makeTableRaw(input_tree,h_counts,cut_string,'rawPFMET/fauxMET > 200',n_tot_entries,counter)

# Jet eta cut
cut_string = cut_string + '&&(abs(jet1.eta()) < 2.5)'
counter = makeTableRaw(input_tree,h_counts,cut_string,'abs(AK5_{eta}) < 2.5',n_tot_entries,counter)

# Njets cut
cut_string = cut_string + '&&(njets < 3)'
counter = makeTableRaw(input_tree,h_counts,cut_string,'nAK5 < 3',n_tot_entries,counter)

## Njets cut
#cut_string = cut_string + '&&(njets == 1)'
#counter = makeTableRaw(input_tree,h_counts,cut_string,'nAK5 == 1',n_tot_entries,counter)

# Delta phi jet1,jet2
cut_string = cut_string + '&&(njets == 1 || TMath::ACos(TMath::Cos(jet2.phi() - jet1.phi())) < 2.0)'
counter = makeTableRaw(input_tree,h_counts,cut_string,'dPhi(AK5_1,AK5_2) < 2.0',n_tot_entries,counter)

# Lepton veto
cut_string = cut_string + '&&(nlep == 0)'
counter = makeTableRaw(input_tree,h_counts,cut_string,'nlep == 0',n_tot_entries,counter)

# Photon veto
if (options.input_selection == 'signal'):
  cut_string = cut_string + '&&(nphotons == 0)'
  counter = makeTableRaw(input_tree,h_counts,cut_string,'nphotons == 0',n_tot_entries,counter)

# Tau veto
cut_string = cut_string + '&&(ntaus == 0)'
counter = makeTableRaw(input_tree,h_counts,cut_string,'ntaus == 0',n_tot_entries,counter)

# Met noise filters
cut_string = cut_string + '&&(metFiltersWord == 511 || metFiltersWord == 1023)'
counter = makeTableRaw(input_tree,h_counts,cut_string,'MET noise filters',n_tot_entries,counter)

# Jet noise filters
cut_string = cut_string + '&&(jet1CHF > 0.2 && jet1NHF < 0.7 && jet1NEMF < 0.7)'
counter = makeTableRaw(input_tree,h_counts,cut_string,'AK5_1 noise filters',n_tot_entries,counter)

print '</verbatim>'

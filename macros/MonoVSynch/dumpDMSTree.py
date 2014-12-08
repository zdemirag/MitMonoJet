#! /usr/bin/env python

import os, re, array, ROOT
import sys
from optparse import OptionParser

# - S U B R O U T I N E S -----------------------------------------------------------------------------

# Print Event information
def printEventInfo(tree):
  print '{0:10s}: {1:10d} | {2:10s}: {3:10d} | {4:10s}: {5:10d}'\
  .format('Run',tree.run,'LumiSec',tree.lumi,'EventNum',tree.event) 
  return

# Print MET information
def printMetInfo(tree):
  print '{0:10s}: {1:10.4f} | {2:10s}: {3:10.4f}'\
  .format('PFMet',tree.metRaw,'MVAMet',tree.mvamet) 
  return

# Print FatJet information
def printFatJetInfo(tree):
  # First Jet
  print '{0:10s}: {1:10.4f} | {2:10s}: {3:10.4f} | {4:10s}: {5:10.4f}'\
  .format('FJet1Pt',tree.fjet1.Pt(),'FJet1Eta',tree.fjet1.Eta(),'FJet1Btag',tree.fjet1Btag) 

  # Second Jet
  print '{0:10s}: {1:10.4f} | {2:10s}: {3:10.4f} | {4:10s}: {5:10.4f}'\
  .format('FJet2Pt',tree.fjet2.Pt(),'FJet2Eta',tree.fjet2.Eta(),'FJet2Btag',tree.fjet2Btag) 
  return

# Print Jet information
def printJetInfo(tree):
  # First Jet
  print '{0:10s}: {1:10.4f} | {2:10s}: {3:10.4f}'\
  .format('Jet1Pt',tree.jet1.Pt(),'Jet1Eta',tree.jet1.Eta()) 

  # Second Jet
  print '{0:10s}: {1:10.4f} | {2:10s}: {3:10.4f}'\
  .format('Jet2Pt',tree.jet2.Pt(),'Jet2Eta',tree.jet2.Eta()) 

  # Third Jet
  print '{0:10s}: {1:10.4f} | {2:10s}: {3:10.4f}'\
  .format('Jet3Pt',tree.jet3.Pt(),'Jet3Eta',tree.jet3.Eta()) 
  return

# Print BJet information
def printBJetInfo(tree):
  # First BJet
  print '{0:10s}: {1:10.4f} | {2:10s}: {3:10.4f} | {4:10s}: {5:10.4f}'\
  .format('BJet1Pt',tree.bjet1.Pt(),'BJet1Eta',tree.bjet1.Eta(),'BJet1Btag',tree.bjet1Btag) 

  # Second BJet
  print '{0:10s}: {1:10.4f} | {2:10s}: {3:10.4f} | {4:10s}: {5:10.4f}'\
  .format('BJet2Pt',tree.bjet2.Pt(),'BJet2Eta',tree.bjet2.Eta(),'BJet2Btag',tree.bjet2Btag) 
  return

# Print Lepton information
def printLepInfo(tree):
  # First Lepton
  print '{0:10s}: {1:10.4f} | {2:10s}: {3:10.4f} | {4:10s}: {5:10d}'\
  .format('Lep1Pt',tree.lep1.Pt(),'Lep1Eta',tree.lep1.Eta(),'Lep1Id',tree.lid1) 

  # Second Lepton
  print '{0:10s}: {1:10.4f} | {2:10s}: {3:10.4f} | {4:10s}: {5:10d}'\
  .format('Lep1Pt',tree.lep2.Pt(),'Lep1Eta',tree.lep2.Eta(),'Lep1Id',tree.lid2) 
  return

# Print Sub information
def printSubInfo(tree):
  # Send a message to the user
  print 'All the following quantities are obtained from the hardest fat jet in the event' 

  # Masses
  print '{0:10s}: {1:10.4f} | {2:10s}: {3:10.4f} | {4:10s}: {5:10.4f} | {6:10s}: {7:10.4f}'\
  .format('Mass',tree.fjet1.M(),'MassPrun',tree.fjet1MassPruned,'MassFilt',tree.fjet1MassFiltered,'MassTrim',tree.fjet1MassTrimmed) 
  print '{0:10s}: {1:10.4f} | {2:10s}: {3:10.4f} | {4:10s}: {5:10.4f} | {6:10s}: {7:10.4f}'\
  .format('MassSDbm1',tree.fjet1MassSDbm1,'MassSDb0',tree.fjet1MassSDb0,'MassSDb1',tree.fjet1MassSDb1,'MassSDb2',tree.fjet1MassSDb2) 

  # Sub-jettiness
  print '{0:10s}: {1:10.4f} | {2:10s}: {3:10.4f} | {4:10s}: {5:10.4f}'\
  .format('Tau1',tree.fjet1Tau1,'Tau2',tree.fjet1Tau2,'Tau3',tree.fjet1Tau3) 

  # ECF
  print '{0:10s}: {1:10.4f} | {2:10s}: {3:10.4f} | {4:10s}: {5:10.4f} | {6:10s}: {7:10.4f}'\
  .format('C2b0',tree.fjet1C2b0,'C2b0p5',tree.fjet1C2b0p5,'C2b1',tree.fjet1C2b1,'C2b2',tree.fjet1C2b2) 

  # Others
  print '{0:10s}: {1:10.4f} | {2:10s}: {3:10.4f} | {4:10s}: {5:10.4f} | {6:10s}: {7:10.4f}'\
  .format('QJetVol',tree.fjet1QJetVol,'QGTag',tree.fjet1QGtag,'Pull',tree.fjet1Pull,'PullAngle',tree.fjet1PullAngle) 

  # Subjets
  print '{0:10s}: {1:10d} | {2:10s}: {3:10.4f} | {4:10s}: {5:10.4f}'\
  .format('NSubJets',tree.fjet1nsj,'QGTagSub1',tree.fjet1QGtagSub1,'QGTagSub2',tree.fjet1QGtagSub2) 
  print '{0:10s}: {1:10.4f} | {2:10s}: {3:10.4f} | {4:10s}: {5:10.4f} | {6:10s}: {7:10.4f}'\
  .format('PtSub1',tree.fjet1sj1.Pt(),'EtaSub1',tree.fjet1sj1.Eta(),'PtSub2',tree.fjet1sj2.Pt(),'EtaSub2',tree.fjet1sj2.Eta()) 

  return

# - M A I N ----------------------------------------------------------------------------------------

# Prepare the command line parser
parser = OptionParser()
parser.add_option("-f", "--file", dest="input_file", default='boostedv_s12-sch-v7a_noskim_flatntuple.root',
                  help="input root file [default: %default]")
parser.add_option("-t", "--treename", dest="input_tree", default='DMSTree',
                  help="root tree name [default: %default]")
parser.add_option("-n", "--nprocs", dest="nprocs", type="int", default=10,
                  help="number of processed entries [default: %default]")
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

# Check the number of entries in the tree
n_entries = input_tree.GetEntriesFast()
print 'INFO - Input tree entries: ' + str(n_entries)

# Determine number of entries to process
if options.nprocs > 0 and options.nprocs < n_entries:
  n_entries = options.nprocs

print 'INFO - Number of entries to be processed: ' + str(n_entries)

# Loop over the entries
for ientry in range(0,n_entries):
  # Grab the n'th entry
  input_tree.GetEntry(ientry)

  print 'INFO ------------------------ New Event ------------------------ '

  # Print event information
  print '\n'
  printEventInfo(input_tree)

  # Print MET information
  print '\n'
  printMetInfo(input_tree)

  # Print FatJet information
  print '\n'
  printFatJetInfo(input_tree)

  # Print Jet information
  print '\n'
  printJetInfo(input_tree)

  # Print BJet information
  print '\n'
  printBJetInfo(input_tree)

  # Print Lepton information
  print '\n'
  printLepInfo(input_tree)

  # Print Substructure information
  print '\n'
  printSubInfo(input_tree)

  print '\n'


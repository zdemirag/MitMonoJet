#! /usr/bin/env python

import os, re, array, ROOT
import sys
from optparse import OptionParser

# - S U B R O U T I N E S -----------------------------------------------------------------------------

# Make table header
def makeTableHeader(region):
  if region == 'Met':
    print ' {0:>9s} & {1:>8s} & {2:>8s} & {3:>7s} & {4:>7s} & {5:>9s} & {6:>7s} & {7:>7s} \\\ \hline'\
    .format('','Z+jets','W+jets','Diboson','Top','Total','VH','ggH') 
  if region == 'Zll':
    print ' {0:>9s} & {1:>7s} & {2:>7s} & {3:>7s} & {4:>7s} & {5:>7s} \\\ \hline'\
    .format('','Z+jets','Diboson','Top','Total','Data') 
  if region == 'Wlv':
    print ' {0:>9s} & {1:>9s} & {2:>7s} & {3:>7s} & {4:>7s} & {5:>9s} & {6:>9s} \\\ \hline'\
    .format('','W+jets','Top','Diboson','Z+jets','Total','Data') 
  if region == 'Pj':
    print ' {0:>9s} & {1:>9s} & {2:>9s} & {3:>9s} & {4:>9s} \\\ \hline'\
    .format('','G+jets','QCD','Total','Data') 

  return

# Make table raw
def makeTableRaw(region,rootfile,rawName=''):
  if region == 'Met':
    # Get the entries
    nZjets = rootfile.Get('Z+jets').GetSumOfWeights()
    nWjets = rootfile.Get('W+jets').GetSumOfWeights()
    nDiboson = rootfile.Get('Diboson').GetSumOfWeights()
    nTop = rootfile.Get('Top').GetSumOfWeights()
    nOthers = rootfile.Get('Others').GetSumOfWeights()
    nWH = rootfile.Get('WH_125').GetSumOfWeights()
    nZH = rootfile.Get('ZH_125').GetSumOfWeights()
    nggH = rootfile.Get('ggH_125').GetSumOfWeights()
    # Prepare the relevant groups
    nTotal = nZjets + nWjets + nDiboson + nTop + nOthers 
    nVH = nWH + nZH
    print ' {0:>9s} & {1:8.2f} & {2:8.2f} & {3:7.2f} & {4:7.2f} & {5:9.2f} & {6:7.2f} & {7:7.2f} \\\ \hline'\
    .format(rawName,nZjets,nWjets,nDiboson,nTop,nTotal,nVH,nggH) 

  if region == 'Zll':
    # Get the entries
    nZjets = rootfile.Get('Z+jets').GetSumOfWeights()
    nWjets = rootfile.Get('W+jets').GetSumOfWeights()
    nDiboson = rootfile.Get('Diboson').GetSumOfWeights()
    nTop = rootfile.Get('Top').GetSumOfWeights()
    nOthers = rootfile.Get('Others').GetSumOfWeights()
    nData = rootfile.Get('htempdata_0').GetSumOfWeights()
    # Prepare the relevant groups
    nTotal = nZjets + nWjets + nDiboson + nTop + nOthers 
    print ' {0:>9s} & {1:7.2f} & {2:7.2f} & {3:7.2f} & {4:7.2f} & {5:7.2f} \\\ \hline'\
    .format(rawName,nZjets,nDiboson,nTop,nTotal,nData) 

  if region == 'Wlv':
    # Get the entries
    nZjets = rootfile.Get('Z+jets').GetSumOfWeights()
    nWjets = rootfile.Get('W+jets').GetSumOfWeights()
    nDiboson = rootfile.Get('Diboson').GetSumOfWeights()
    nTop = rootfile.Get('Top').GetSumOfWeights()
    nOthers = rootfile.Get('Others').GetSumOfWeights()
    nData = rootfile.Get('htempdata_0').GetSumOfWeights()
    # Prepare the relevant groups
    nTotal = nZjets + nWjets + nDiboson + nTop + nOthers 
    print ' {0:>9s} & {1:9.2f} & {2:7.2f} & {3:7.2f} & {4:7.2f} & {5:9.2f} & {6:9.2f} \\\ \hline'\
    .format(rawName,nWjets,nTop,nDiboson,nZjets,nTotal,nData) 

  if region == 'Pj':
    # Get the entries
    nGjets = rootfile.Get('G+jets').GetSumOfWeights()
    nZjets = rootfile.Get('Z+jets').GetSumOfWeights()
    nWjets = rootfile.Get('W+jets').GetSumOfWeights()
    nDiboson = rootfile.Get('Diboson').GetSumOfWeights()
    nTop = rootfile.Get('Top').GetSumOfWeights()
    nOthers = rootfile.Get('Others').GetSumOfWeights()
    nData = rootfile.Get('htempdata_0').GetSumOfWeights()
    # Prepare the relevant groups
    nTotal = nGjets +  nZjets + nWjets + nDiboson + nTop + nOthers 
    print ' {0:>9s} & {1:9.2f} & {2:9.2f} & {3:9.2f} & {4:9.2f} \\\ \hline'\
    .format(rawName,nGjets,nOthers,nTotal,nData) 

  return
    
# - M A I N ----------------------------------------------------------------------------------------

# Prepare the command line parser
parser = OptionParser()
parser.add_option("-f", "--file", dest="input_folder", default='../MonoVPlots/',
                  help="input folder [default: %default]")
parser.add_option("-r", "--region", dest="region", default='Met',
                  help="region to analyze [default: %default]")
(options, args) = parser.parse_args()

# Get all the root classes
from ROOT import *

# Signal region table
if options.region == 'Met':
  print 'INFO - Creating signal region tables.'
  makeTableHeader(options.region)

  input_file = TFile.Open(options.input_folder + 'BDT_Met_baseline_isData.root')
  if not input_file:
    print 'ERROR - Cannot open input root file: ' + options.input_file + ' , exiting!'
    raise SystemExit
  makeTableRaw(options.region,input_file,'baseline')
  input_file = TFile.Open(options.input_folder + 'BDT_Met_inclusive_isData.root')
  if not input_file:
    print 'ERROR - Cannot open input root file: ' + options.input_file + ' , exiting!'
    raise SystemExit
  makeTableRaw(options.region,input_file,'inclusive')
  input_file = TFile.Open(options.input_folder + 'BDT_Met_isData.root')
  if not input_file:
    print 'ERROR - Cannot open input root file: ' + options.input_file + ' , exiting!'
    raise SystemExit
  makeTableRaw(options.region,input_file,'boosted')

# Zll region table
if options.region == 'Zll':
  print 'INFO - Creating Zll control region tables.'
  makeTableHeader(options.region)

  input_file = TFile.Open(options.input_folder + 'BDT_Zll_baseline_isData.root')
  if not input_file:
    print 'ERROR - Cannot open input root file: ' + options.input_file + ' , exiting!'
    raise SystemExit
  makeTableRaw(options.region,input_file,'baseline')
  input_file = TFile.Open(options.input_folder + 'BDT_Zll_inclusive_isData.root')
  if not input_file:
    print 'ERROR - Cannot open input root file: ' + options.input_file + ' , exiting!'
    raise SystemExit
  makeTableRaw(options.region,input_file,'inclusive')
  input_file = TFile.Open(options.input_folder + 'BDT_Zll_isData.root')
  if not input_file:
    print 'ERROR - Cannot open input root file: ' + options.input_file + ' , exiting!'
    raise SystemExit
  makeTableRaw(options.region,input_file,'boosted')

# Wlv region table
if options.region == 'Wlv':
  print 'INFO - Creating Wlv control region tables.'
  makeTableHeader(options.region)

  input_file = TFile.Open(options.input_folder + 'BDT_Wlv_baseline_isData.root')
  if not input_file:
    print 'ERROR - Cannot open input root file: ' + options.input_file + ' , exiting!'
    raise SystemExit
  makeTableRaw(options.region,input_file,'baseline')
  input_file = TFile.Open(options.input_folder + 'BDT_Wlv_inclusive_isData.root')
  if not input_file:
    print 'ERROR - Cannot open input root file: ' + options.input_file + ' , exiting!'
    raise SystemExit
  makeTableRaw(options.region,input_file,'inclusive')
  input_file = TFile.Open(options.input_folder + 'BDT_Wlv_isData.root')
  if not input_file:
    print 'ERROR - Cannot open input root file: ' + options.input_file + ' , exiting!'
    raise SystemExit
  makeTableRaw(options.region,input_file,'boosted')

# Wlv region table
if options.region == 'Pj':
  print 'INFO - Creating Photon+jets control region tables.'
  makeTableHeader(options.region)

  input_file = TFile.Open(options.input_folder + 'BDT_Pj_baseline_isData.root')
  if not input_file:
    print 'ERROR - Cannot open input root file: ' + options.input_file + ' , exiting!'
    raise SystemExit
  makeTableRaw(options.region,input_file,'baseline')
  input_file = TFile.Open(options.input_folder + 'BDT_Pj_inclusive_isData.root')
  if not input_file:
    print 'ERROR - Cannot open input root file: ' + options.input_file + ' , exiting!'
    raise SystemExit
  makeTableRaw(options.region,input_file,'inclusive')
  input_file = TFile.Open(options.input_folder + 'BDT_Pj_isData.root')
  if not input_file:
    print 'ERROR - Cannot open input root file: ' + options.input_file + ' , exiting!'
    raise SystemExit
  makeTableRaw(options.region,input_file,'boosted')
  

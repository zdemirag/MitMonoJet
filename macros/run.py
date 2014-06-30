#!/bin/env python
# -*- coding: utf-8 -*-
'''
run.py
======
Wrapper around a run macro to ease passing the arguments.

It replaces the following command in ROOT:

    .x runMonoJet.C+("0000", "noskim", "s12-h125inv-gf-v7a", "filefi/032",
                     "/home/cmsprod/catalog/t2mit", "monojet-2013-032", 1)

By calling this from the command line (after editing the parameters here):

    ./run.py

where run.py contains:

    runMacro   = 'runMonoJet.C'
    skim       = 'noskim'
    catalogDir = '/home/cmsprod/catalog/t2mit'
    outputName = 'monojet-2014-032'

    dataset    = 's12-h125inv-gf-v7a'
    book       = 'filefi/032'
    fileset    = '0000'
    nEvents    = 1

Jan Veverka, 20 June 2014
'''

import sys
sys.argv.append('-b') ## Load ROOT in batch mode
import ROOT

## Monojet Stuff
runMacro   = 'runMonoJet.C'
skim       = 'noskim'
catalogDir = '/home/cmsprod/catalog/t2mit'
outputName = 'photon-id'

dataset    = 's12-pj80_120-v7a'
book       = 'filefi/032'
fileset    = '0000'
nEvents    = 100

ROOT.gROOT.ProcessLine('.L ' + runMacro + '++' )
run = getattr(ROOT, runMacro.split('.')[0])
run(fileset, skim, dataset, book, catalogDir, outputName, nEvents)

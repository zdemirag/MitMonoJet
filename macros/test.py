#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Scan the content of the tree to quickly check that its production worked fine.
USAGE: ./test.py
'''

import glob
import os
import sys
## Set batch mode for root
sys.argv.append('-b')

import ROOT

#_______________________________________________________________________________
def get_latest_root_file(mask = '*.root'):
    filenames =  glob.glob(mask)
    mtimes = [os.stat(f)[-2] for f in filenames]
    return max(zip(mtimes, filenames))[1]
## get_latest_root_file


filename = get_latest_root_file('photon-id_*_noskim_*.root')
dirname = 'MonoJetTreeWriter'
treename = 'MJetTree'


#_______________________________________________________________________________
def main():
    ## Make the tree global so that it can be accessed in an interactive
    ## session.
    global tree
    print '%s:%s/%s' % (filename, dirname, treename)
    tree = get_tree(filename, dirname, treename)
    # print_all_branches(tree)
    tree.Scan('run:lumi:event:pho1.Pt():jet1.Pt()', 'nphotons>0', '', 20)
## main


#_______________________________________________________________________________
def get_tree(filename, dirname, treename):
    ## Make the source file object global to prevent its garbage collection,
    ## which destroys the tree as a side effect.
    global source
    source = ROOT.TFile(filename)
    tree = source.FindObjectAny(dirname).Get(treename)
    return tree
## get_tree


#_______________________________________________________________________________
def print_all_branches(tree, separator = '   ', width = 80):
    names = [b.GetName() for b in tree.GetListOfBranches()]
    maxl = max([len(n) for n in names])
    sepl = len(separator)
    ncolumns = width / (maxl + sepl)
    nlines = len(names) / ncolumns + 1
    for line in range(nlines):
        entries = []
        for column in range(ncolumns):
            index = line + nlines * column
            if index < len(names):
                entries.append(names[index].ljust(maxl))
        print separator.join(entries)
## print_all_branches


#______________________________________________________________________________
if __name__ == '__main__':
    import user
    main()

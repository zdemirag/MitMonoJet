#!/usr/bin/env python

from ROOT import TFile

import argparse
parser = argparse.ArgumentParser(description='adjusting crazy bambu structure')
parser.add_argument('-i', action='store', dest='input', help='input file')
parser.add_argument('-o', action='store', dest='output', help='dataset')
args = parser.parse_args()

if not args.input: args.input='/home/mzanetti/cms/hist/monojet-2013-June12/filefi/031/s12-h125inv-vbf-v7a/monojet-2013-June12_s12-h125inv-vbf-v7a_noskim_0000.root'
input = TFile(args.input)
hDAllEvents = input.FindObjectAny('hDAllEvents')
hNPUTrue = input.FindObjectAny('hNPUTrue')

if not args.output: args.output='tmp.root'
output = TFile(args.output,'recreate')
hDAllEvents.Write()
hNPUTrue.Write()

tree_names = {'MJetTree':'tree','MJetTreeDiLepton':'treeZll','MJetTreeWlnu':'treeWln',}
trees, newTrees = {}, {} 
for tree_name in tree_names:
    input.cd()
    trees[tree_name] = input.FindObjectAny(tree_name)
    trees[tree_name].SetName(tree_names[tree_name])
    output.cd()
    newTrees[tree_name] = trees[tree_name].CopyTree('')
    newTrees[tree_name].Write()
#output.Write()

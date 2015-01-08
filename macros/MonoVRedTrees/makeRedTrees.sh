#!/bin/bash
#---------------------------------------------------------------------------------------------------
# Make the reduced trees for limit computation for all regions
#
#                                                                         L.Di Matteo (Dec 14, 2014)
#---------------------------------------------------------------------------------------------------
# Boosted category
root -b -q ../rootlogon_monojet.C makeReducedTree.C+'(0,19700,false)'
root -b -q ../rootlogon_monojet.C makeReducedTree.C+'(1,19700,true)'
root -b -q ../rootlogon_monojet.C makeReducedTree.C+'(2,19700,true)'
root -b -q ../rootlogon_monojet.C makeReducedTree.C+'(3,19700,true)'
# Resolved category
root -b -q ../rootlogon_monojet.C makeReducedTree.C+'(4,19700,false)'
root -b -q ../rootlogon_monojet.C makeReducedTree.C+'(5,19700,true)'
root -b -q ../rootlogon_monojet.C makeReducedTree.C+'(6,19700,true)'
root -b -q ../rootlogon_monojet.C makeReducedTree.C+'(7,19700,true)'
# Inclusive category
root -b -q ../rootlogon_monojet.C makeReducedTree.C+'(8,19700,false)'
root -b -q ../rootlogon_monojet.C makeReducedTree.C+'(9,19700,true)'
root -b -q ../rootlogon_monojet.C makeReducedTree.C+'(10,19700,true)'
root -b -q ../rootlogon_monojet.C makeReducedTree.C+'(11,19700,true)'
# Baseline selection (for synch), i.e. no exclusivity in the selection
#root -b -q ../rootlogon_monojet.C makeReducedTree.C+'(8,19700,false,false)'
#root -b -q ../rootlogon_monojet.C makeReducedTree.C+'(9,19700,true,false)'
#root -b -q ../rootlogon_monojet.C makeReducedTree.C+'(10,19700,true,false)'
#root -b -q ../rootlogon_monojet.C makeReducedTree.C+'(11,19700,true,false)'

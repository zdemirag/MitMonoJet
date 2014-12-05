#!/bin/bash
#---------------------------------------------------------------------------------------------------
# Make the reduced trees for limit computation for all regions
#
#                                                                         L.Di Matteo (Dec 14, 2014)
#---------------------------------------------------------------------------------------------------
root -b -q ../rootlogon_monojet.C makeReducedTree.C+'(0,19700,false)'
root -b -q ../rootlogon_monojet.C makeReducedTree.C+'(1,19700,true)'
root -b -q ../rootlogon_monojet.C makeReducedTree.C+'(2,19700,true)'
root -b -q ../rootlogon_monojet.C makeReducedTree.C+'(3,19700,true)'


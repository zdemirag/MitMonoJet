#!/bin/bash
#---------------------------------------------------------------------------------------------------
# Prepare weighted histogram for nice plotting
#
#                                                                         L.Di Matteo (Nov 14, 2014)
#---------------------------------------------------------------------------------------------------
# signal region plots
root -b -q ../rootlogon_monojet.C makeRawPlot.C+'(19700,1,"metRaw",4,55,250,800)'
root -b -q ../rootlogon_monojet.C makeRawPlot.C+'(19700,1,"fjet1.M()",4,40,0,200)'
root -b -q ../rootlogon_monojet.C makeRawPlot.C+'(19700,1,"fjet1.Pt()",4,55,250,800)'
root -b -q ../rootlogon_monojet.C makeRawPlot.C+'(19700,1,"bdt_all",5,40,-1.,1.)'

exit 0

#!/bin/bash
#---------------------------------------------------------------------------------------------------
# Make the reduced trees for limit computation for all regions
#
#                                                                         L.Di Matteo (Dec 14, 2014)
#---------------------------------------------------------------------------------------------------
if [ $1 = "boosted" ]  || [ $1 = "all" ];
then
  # Boosted category
  root -b -q ../rootlogon_monojet.C makeReducedTree.C+'(0,19700,false,true,false,0)'
  root -b -q ../rootlogon_monojet.C makeReducedTree.C+'(1,19700,true,true,false,0)'
  root -b -q ../rootlogon_monojet.C makeReducedTree.C+'(2,19700,true,true,false,0)'
  root -b -q ../rootlogon_monojet.C makeReducedTree.C+'(3,19700,true,true,false,0)'
fi
if [ $1 = "resolved" ] || [ $1 = "all" ];
then
  # Resolved category
  root -b -q ../rootlogon_monojet.C makeReducedTree.C+'(4,19700,false,true,false,0)'
  root -b -q ../rootlogon_monojet.C makeReducedTree.C+'(5,19700,true,true,false,0)'
  root -b -q ../rootlogon_monojet.C makeReducedTree.C+'(6,19700,true,true,false,0)'
  root -b -q ../rootlogon_monojet.C makeReducedTree.C+'(7,19700,true,true,false,0)'
fi
if [ $1 = "inclusive" ] || [ $1 = "all" ];
  then
  # Inclusive category
  root -b -q ../rootlogon_monojet.C makeReducedTree.C+'(8,19700,false,true,false,0)'
  root -b -q ../rootlogon_monojet.C makeReducedTree.C+'(9,19700,true,true,false,0)'
  root -b -q ../rootlogon_monojet.C makeReducedTree.C+'(10,19700,true,true,false,0)'
  root -b -q ../rootlogon_monojet.C makeReducedTree.C+'(11,19700,true,true,false,0)'
fi

#JES Up
if [ $1 = "boostedJESup" ]  || [ $1 = "allJESup" ];
then
  # Boosted category
  root -b -q ../rootlogon_monojet.C makeReducedTree.C+'(0,19700,false,true,false,1)'
  root -b -q ../rootlogon_monojet.C makeReducedTree.C+'(1,19700,true,true,false,1)'
  root -b -q ../rootlogon_monojet.C makeReducedTree.C+'(2,19700,true,true,false,1)'
  root -b -q ../rootlogon_monojet.C makeReducedTree.C+'(3,19700,true,true,false,1)'
fi
if [ $1 = "resolvedJESup" ] || [ $1 = "allJESup" ];
then
  # Resolved category
  root -b -q ../rootlogon_monojet.C makeReducedTree.C+'(4,19700,false,true,false,1)'
  root -b -q ../rootlogon_monojet.C makeReducedTree.C+'(5,19700,true,true,false,1)'
  root -b -q ../rootlogon_monojet.C makeReducedTree.C+'(6,19700,true,true,false,1)'
  root -b -q ../rootlogon_monojet.C makeReducedTree.C+'(7,19700,true,true,false,1)'
fi
if [ $1 = "inclusiveJESup" ] || [ $1 = "allJESup" ];
  then
  # Inclusive category
  root -b -q ../rootlogon_monojet.C makeReducedTree.C+'(8,19700,false,true,false,1)'
  root -b -q ../rootlogon_monojet.C makeReducedTree.C+'(9,19700,true,true,false,1)'
  root -b -q ../rootlogon_monojet.C makeReducedTree.C+'(10,19700,true,true,false,1)'
  root -b -q ../rootlogon_monojet.C makeReducedTree.C+'(11,19700,true,true,false,1)'
fi

#JES Down
if [ $1 = "boostedJESdown" ]  || [ $1 = "allJESdown" ];
then
  # Boosted category
  root -b -q ../rootlogon_monojet.C makeReducedTree.C+'(0,19700,false,true,false,-1)'
  root -b -q ../rootlogon_monojet.C makeReducedTree.C+'(1,19700,true,true,false,-1)'
  root -b -q ../rootlogon_monojet.C makeReducedTree.C+'(2,19700,true,true,false,-1)'
  root -b -q ../rootlogon_monojet.C makeReducedTree.C+'(3,19700,true,true,false,-1)'
fi
if [ $1 = "resolvedJESdown" ] || [ $1 = "allJESdown" ];
then
  # Resolved category
  root -b -q ../rootlogon_monojet.C makeReducedTree.C+'(4,19700,false,true,false,-1)'
  root -b -q ../rootlogon_monojet.C makeReducedTree.C+'(5,19700,true,true,false,-1)'
  root -b -q ../rootlogon_monojet.C makeReducedTree.C+'(6,19700,true,true,false,-1)'
  root -b -q ../rootlogon_monojet.C makeReducedTree.C+'(7,19700,true,true,false,-1)'
fi
if [ $1 = "inclusiveJESdown" ] || [ $1 = "allJESdown" ];
  then
  # Inclusive category
  root -b -q ../rootlogon_monojet.C makeReducedTree.C+'(8,19700,false,true,false,-1)'
  root -b -q ../rootlogon_monojet.C makeReducedTree.C+'(9,19700,true,true,false,-1)'
  root -b -q ../rootlogon_monojet.C makeReducedTree.C+'(10,19700,true,true,false,-1)'
  root -b -q ../rootlogon_monojet.C makeReducedTree.C+'(11,19700,true,true,false,-1)'
fi

# Only for testing
if [ $1 = "test" ];
  then
  root -b -q ../rootlogon_monojet.C makeReducedTree.C+'(10,19700,false,true,true,0)'
fi
if [ $1 = "testJESup" ];
  then
  root -b -q ../rootlogon_monojet.C makeReducedTree.C+'(10,19700,false,true,true,1)'
fi
if [ $1 = "testJESdown" ];
  then
  root -b -q ../rootlogon_monojet.C makeReducedTree.C+'(10,19700,false,true,true,-1)'
fi

exit 0

#!/bin/bash
#---------------------------------------------------------------------------------------------------
# Make the tables for all the regions
#
#                                                                         L.Di Matteo (Dec 20, 2014)
#---------------------------------------------------------------------------------------------------
./makeMonoVTables.py -r Met >  tables.tex
./makeMonoVTables.py -r Zll >> tables.tex
./makeMonoVTables.py -r Wlv >> tables.tex
./makeMonoVTables.py -r Pj  >> tables.tex

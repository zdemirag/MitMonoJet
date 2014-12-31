#! /bin/bash
#===================================================================================================
#
# Script to add the MVA bdt branch to the flattrees. It can be used as model for other additions.
# The only input parameter, $1$, should be unset only if the entire flattrees have to be patched.
#
#                                                                         L.Di Matteo (Nov 11, 2014)
#===================================================================================================

mkdir -p /tmp/$USER

echo $1

if [ -z $1 ]
then 
  for file in `ls /mnt/hscratch/$USER/$MIT_PROD_CFG/merged/*_flatntuple.root`; do
    echo "INFO - Processing: $file"
    mv $file /tmp/$USER/in.root
    root -q -b computeBDT.C+\(\"/tmp/$USER/in.root\"\)
    root -q -b patchTree.C+\(\"/tmp/$USER/in.root\"\)
    mv /tmp/$USER/in.root $file
  done
else 
  for dataset in `cat ${MIT_USER_DIR}/config/${MIT_PROD_CFG}.txt|grep -v ^#|tr -s ' '|cut -d' ' -f 2`; do
    file="/mnt/hscratch/$USER/${MIT_PROD_CFG}/merged/${MIT_PROD_CFG}_${dataset}_noskim_flatntuple.root"
    echo "INFO - Processing: $file"
    mv $file /tmp/$USER/in.root
    root -q -b computeBDT.C+\(\"/tmp/$USER/in.root\"\)
    root -q -b patchTree.C+\(\"/tmp/$USER/in.root\"\)
    mv /tmp/$USER/in.root $file
  done
fi



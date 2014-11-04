#! /bin/bash

mkdir -p /tmp/$USER

for file in `ls /mnt/hscratch/$USER/$MIT_PROD_CFG/merged/*_flatntuple.root`; do
    echo "INFO - Processing: $file"
    mv $file /tmp/$USER/in.root
    root -q -b computeBDT.C+\(\"/tmp/$USER/in.root\"\)
    root -q -b patchTree.C+\(\"/tmp/$USER/in.root\"\)
    mv /tmp/$USER/in.root $file
done


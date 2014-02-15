#!/bin/bash
#===================================================================================================
# Script to merge all histogram files from a complete analysis task such that each sample with
# a well defined cross section is in one 'merged' file.
#
#                                                                             Ch.Paus (Aug 15, 2010)
#===================================================================================================
PATTERN="$1"

echo " Config: ${MIT_MONOJET_DIR}/config/${MIT_PROD_CFG}.txt"

# make clean copy of the config
grep -v ^# ${MIT_MONOJET_DIR}/config/${MIT_PROD_CFG}.txt > /tmp/mergeMonoJet.$$
nLine=0
rm processingStatus.pkl

for dset in `cat /tmp/mergeMonoJet.$$|grep -v ^#|tr -s ' '|cut -d' ' -f 2`
do
  nLine=$(($nLine+1))

  # find the line to this dataset and do further analysis
  line=`sed -n ${nLine}p /tmp/mergeMonoJet.$$`

  # determine the input dataset
  BOOK_VERSION=`echo $line | tr -s ' ' | cut -d ' ' -f 1`
  DATASET=`     echo $line | tr -s ' ' | cut -d ' ' -f 2`
  SKIM=`        echo $line | tr -s ' ' | cut -d ' ' -f 3`

  if [ -z $DATASET ]; then continue; fi

  # adjust book / catalog dir
  BOOK=`echo $BOOK_VERSION | cut -d/ -f2-3`
  CEXT=`echo $BOOK_VERSION | cut -d/ -f1`


  # skip merging if there are still some jobs running
  job_counter=`condor_q mzanetti -w | grep -c $DATASET`; 
  if [ $job_counter -gt 1 ] ; then
      echo "$job_counter jobs still running for $DATASET. STILL RUNNING"
      ./buildNtupleStatus.py RUNNING $DATASET
      continue
  else
      
      # check if files are there
      if [ ! -d $MIT_PROD_HIST/$MIT_PROD_CFG/$BOOK/${DATASET} ]; then
	  #echo "neither jobs running nor ouput for ${DATASET}. PROCESSING NOT STARTED YET"
	  ./buildNtupleStatus.py NOT_PROCESSED $DATASET
	  continue
      else
	  # check if merged file is already there and ask what to do with the dataset
	  if [ -e $MIT_PROD_HIST/$MIT_PROD_CFG/merged/${DATASET}.root ]; then
	      #echo "merged file for ${DATASET} exists already. DONE"
	      ./buildNtupleStatus.py DONE $DATASET	      
	  else
	      #echo "${DATASET} ready to be merged. TO BE MERGED"
	      ./buildNtupleStatus.py MERGEABLE $DATASET
	  fi
      fi
  fi

done

./printNtupleStatus.py

rm /tmp/mergeMonoJet.$$

exit 0

#!/bin/bash
#===================================================================================================
# Script to merge all histogram files from a complete analysis task such that each sample with
# a well defined cross section is in one 'merged' file.
#
#                                                                             Ch.Paus (Aug 15, 2010)
#===================================================================================================
LABEL="$1"
PATTERN="$2"

echo " Config: ${MIT_MONOJET_DIR}/config/${MIT_PROD_CFG}.txt"

# make clean copy of the config
grep -v ^# ${MIT_MONOJET_DIR}/config/${MIT_PROD_CFG}.txt > /tmp/mergeMonoJet.$$
nLine=0

for dset in `cat /tmp/mergeMonoJet.$$|grep -v ^#|tr -s ' '|cut -d' ' -f 2`
do

  nLine=$(($nLine+1))

  # find the line to this dataset and do further analysis
  line=`sed -n ${nLine}p /tmp/mergeMonoJet.$$`

  # determine the input dataset
  BOOK_VERSION=`echo $line | tr -s ' ' | cut -d ' ' -f 1`
  DATASET=`     echo $line | tr -s ' ' | cut -d ' ' -f 2`
  SKIM=`        echo $line | tr -s ' ' | cut -d ' ' -f 3`

  # adjust book / catalog dir
  BOOK=`echo $BOOK_VERSION | cut -d/ -f2-3`
  CEXT=`echo $BOOK_VERSION | cut -d/ -f1`

  # input and output dir
  INPUT_DIR=$MIT_PROD_HIST/$MIT_PROD_CFG/$BOOK/$LABEL/
  OUTPUT_DIR=$MIT_PROD_HIST/$MIT_PROD_CFG/merged/$LABEL/
  if [ ! -d "$OUTPUT_DIR" ]; then mkdir $OUTPUT_DIR; fi

  # skip merging if there are still some jobs running
  job_counter=`condor_q mzanetti -w | grep -c $DATASET`; 
  #job_counter=0 # FIXME!!
  if [ $job_counter -gt 1 ] ; then
      echo "$job_counter jobs still running for $DATASET, quitting"
      continue

  else

      if [ "$PATTERN" == "" ] || [ "`echo $DATASET | grep $PATTERN`" != "" ]
	  then

          # check if files are there
	  echo $INPUT_DIR/$DATASET
	  if [ ! -d $INPUT_DIR/$DATASET ]; then
	      echo "neither jobs running nor ouput for ${DATASET}"
	      continue
	  fi
	  
          # check if merged file is already there and ask what to do with the dataset
	  #if [ -e $OUTPUT_DIR/${MIT_PROD_CFG}_${DATASET}_${SKIM}.root ]; then
	  if ls $OUTPUT_DIR/${MIT_PROD_CFG}_${DATASET}_${SKIM}*.root &> /dev/null; then
	  #if [ -e $OUTPUT_DIR/${DATASET}.root ]; then
	      read -p "$OUTPUT_DIR/${DATASET}.root exist, should it be deleted? [Y,n]" action
	      if [ "$action" == "n" ]; then
		  echo "ok, skipping";
		  continue
	      fi
	  fi
      
          # now merge the sucker
	  mergeHistSmart.py --Dataset=$DATASET --Skim=$SKIM --InputPath=$INPUT_DIR/ --OutputPath=$OUTPUT_DIR/ --FilenameHeader=$MIT_PROD_CFG

	  # clean up the data (MC done later on)
	  if [ "`echo $DATASET | grep r12`" != "" ] # WARNING!!: "_0" hardcoded in the file name. If more that 1 datafile is produced, the next ones will be ignored
	      then
	      echo 'input for cleaning' $OUTPUT_DIR/${MIT_PROD_CFG}_${DATASET}_${SKIM}_0.root -o $OUTPUT_DIR/${DATASET}_0.root
	      ./cleanMergedOutput.py -i $OUTPUT_DIR/${MIT_PROD_CFG}_${DATASET}_${SKIM}_0.root -o $OUTPUT_DIR/${DATASET}_0.root
	  fi

      fi
  fi

done

rm /tmp/mergeMonoJet.$$

exit 0

#!/usr/bin/python
#---------------------------------------------------------------------------------------------------
# This script can be used to submit condor jobs running on Bavanti nutples to produce flat trees
# this will become obsolete when Bavanti will be included into cmsprod
#---------------------------------------------------------------------------------------------------

import sys, os, shutil

# - S U B R O U T I N E S -----------------------------------------------------------------------------
def makeCondorFile(outputFile,runFile,logDir,datasetName,outRootFileName):
  # Prepare condor submit file  

  line =  'Environment = \"HOSTNAME=' + os.environ['HOSTNAME'] + '\"'
  print >> outputFile, line

  block_str = """\
Universe = vanilla
Requirements = (UidDomain == "mit.edu") && Arch == "X86_64" && Disk >= DiskUsage && (Memory * 1024) >= ImageSize && HasFileTransfer
Notification = Error"""
  # Print the text block into the outputFile file
  print >> outputFile, block_str

  exe_str = 'Executable = ' + runFile
  print >> outputFile, exe_str

  block_str = """\
Rank = Mips
GetEnv = True
Input = /dev/null"""
  # Print the text block into the outputFile file
  print >> outputFile, block_str
  
  line =  'Output = ' + logDir + '/' + datasetName + '.out'
  print >> outputFile, line
  line =  'Error = ' + logDir + '/' + datasetName + '.err'
  print >> outputFile, line
  line =  'Log = ' + logDir + '/' + datasetName + '.log'
  print >> outputFile, line
  line =  'transfer_output_files = ' + outRootFileName
  print >> outputFile, line

  block_str = """\
should_transfer_files = YES
when_to_transfer_output = ON_EXIT"""
  # Print the text block into the outputFile file
  print >> outputFile, block_str
  
  line =  '+AccountingGroup = \"group_cmsuser.' + os.environ['USER'] + '\"'
  print >> outputFile, line
  line =  'Queue'
  print >> outputFile, line
  
  return


# - M A I N -------------------------------------------------------------------------------------------


# => define the number of events you want to process (useful for testing)
nEvents = "-1"

# => turn on testing mode (no submission)
testingMode = False
if testingMode:
  print " WARNING - testing mode on --> jobs will NOT be submitted"

# => save the current directory position
initDir = os.environ['PWD']

# => define the input configuration file using shell environment
inputConfigName = os.environ['MIT_USER_DIR'] + '/config/' + os.environ['MIT_PROD_CFG'] + '.txt' 
print " INFO - preparing flat ntuple submission with configuration file: " + inputConfigName

# => prepare the submit working area
workingDir = os.environ['MIT_PROD_HIST'] + '/' + os.environ['MIT_PROD_CFG'] + '/merged-dev'
if os.path.isdir(workingDir):
  # cleanup all sub and run files in the working area
  os.system('rm ' + workingDir + '/*.sub')
  os.system('rm ' + workingDir + '/*.run')
else:
  os.makedirs(workingDir)
print " INFO - working directory is located under: " + workingDir

# => spell out the output area
outputDir = '/mnt/hscratch/' + os.environ['USER'] + '/' + os.environ['MIT_PROD_CFG'] + '/merged'
print " INFO - output target is: " + outputDir

# => prepare the condor log area
condorLogDir = os.environ['MIT_PROD_LOGS'] + '/' + os.environ['MIT_PROD_CFG'] + '/flat'
if os.path.isdir(condorLogDir):
  shutil.rmtree(condorLogDir)
os.makedirs(condorLogDir)
print " INFO - condor log directory is located under: " + condorLogDir

# => prepare the commands for file transfers
transferFileName = 'flat_ntuple_transfers.txt' 
if os.path.isfile(transferFileName):
  os.system('rm ' + transferFileName)
os.system('touch ' + transferFileName)
transferFile = open(initDir + '/' + transferFileName, 'w')

print " INFO - setup completed, moving now to submissions\n\n "

# move into the working directory
os.chdir(workingDir)
  
# => now read the configuration file line by line
for line in open(inputConfigName):
  li = line.strip()
  # avoid reading commented lines
  if not li.startswith("#"):
    li_array = li.split()      
    # get the relevant parameters for each dataset
    datasetName = li_array[1]
    JSONName = li_array[7]
    if JSONName == "~":
      print " INFO ---> processing MC sample: " + datasetName
    else:
      print " INFO ---> processing dataset: " + datasetName
    # prepare the run instructions for each dataset and put them in a file
    runFile = open(workingDir + '/' + datasetName + ".run", 'w')
    runFile.write('#!/bin/bash\n')
    # first make sure the input file list in created
    runFile.write('ls $MIT_PROD_HIST/$MIT_PROD_CFG/$MIT_PROD_BOOK/032/' + datasetName + '/*ntuple*.root > inputBavanti.txt \n')
    runFile.write('export MIT_PROD_JSON=\'' + JSONName + '\'\n')
    rootString = 'root -b -l -q runFlatBoostedV.C+\\(\\\"000\\\",\\\"noskim\\\",\\\"'
    rootString += datasetName
    rootString += '\\\",\\\"filefi/032\\\",\\\"/home/cmsprod/catalog\\\",\\\"'
    rootString += os.environ['MIT_PROD_CFG']
    rootString += '\\\",' + nEvents + '\\)\n'    
    runFile.write(rootString)
    runFile.close()
    os.chmod(workingDir + '/' + datasetName + ".run", 0755)
    # prepare the condor log dir
    thisCondorLogDir = condorLogDir + '/' + datasetName
    if os.path.isdir(thisCondorLogDir):
      shutil.rmtree(thisCondorLogDir)
    os.makedirs(thisCondorLogDir)
    # prepare the condor submit file
    condorSubFile = open(workingDir + '/' + datasetName + '.sub', 'w')
    thisRootOutFile = os.environ['MIT_PROD_CFG'] + '_' + datasetName + '_noskim_flatntuple.root'
    makeCondorFile(condorSubFile,workingDir + '/' + datasetName + ".run",thisCondorLogDir,datasetName,thisRootOutFile)
    condorSubFile.close()
    # prepare the copy command for final ntuple transfer
    transferFile.write('cp ' + workingDir + '/' + thisRootOutFile + ' ' + outputDir + '\n')
    # submit the jobs if not in testing mode
    if not testingMode:
      os.system('condor_submit ' + workingDir + '/' + datasetName + '.sub')

os.chdir(initDir)

print " \n"
print " INFO - to move the output ntuples to the final destination [" + outputDir + "]:" 
print " INFO - 1) check that the jobs are all finished successfully"
print " INFO - 2) run: source " + transferFileName
transferFile.close()


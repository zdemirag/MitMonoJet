#!/usr/bin/env python

from ROOT import TFile
import sys, os, glob, pickle, argparse
from array import array

parser = argparse.ArgumentParser(description='adjusting crazy bambu structure')
parser.add_argument('-t', action='store', dest='tag', help='tag')
args = parser.parse_args()
if not args.tag: args.tag = ''

# variables
configuration = os.environ['MIT_PROD_CFG']
merged_dir=os.environ['MIT_PROD_HIST']+'/'+os.environ['MIT_PROD_CFG']+'/merged/'+args.tag+'/'
if not os.path.exists(merged_dir): os.makedirs(merged_dir)
output_dir=os.environ['MIT_PROD_HIST']+'/'+os.environ['MIT_PROD_CFG']+'/ntuples/'+args.tag+'/'
if not os.path.exists(output_dir): os.makedirs(output_dir)

print 'Output dir:', output_dir

# scale factors
muonSFs = pickle.load(open('/home/mzanetti/cms/ntupleAnalysis/DM/muonSF.pkl'))

# process config file
xsecs = {}
for line in open(os.environ['MIT_MONOJET_DIR']+'/config/'+configuration+'.txt'):
    # skip commented data
    if line.strip()[0]=='#': continue
    # skip data
    if line.strip().split()[-1]!='~': continue
    xsecs[line.strip().split()[1]] = float(line.strip().split()[4])

# process the datasets
for d in glob.glob(merged_dir+'*.root'):
    dataset = d.replace(merged_dir,'').replace(configuration+'_','').replace('_noskim.root','')
    if not dataset in xsecs.keys(): continue

    print 'processing', dataset 

    ref_file = TFile(d)
    hDAllEvents = ref_file.FindObjectAny('hDAllEvents')
    hNPUTrue = ref_file.FindObjectAny('hNPUTrue')

    # normalization
    ref_initialEntries = hDAllEvents.GetBinContent(1)
    # xsec
    ref_xsec = xsecs[dataset]

    # pileup data
    puFile = TFile('/home/mzanetti/cms/ntupleAnalysis/DM/MyDataPileupHistogram.root')
    pu_target = puFile.FindObjectAny('pileup')
    pu_target.Scale(1.0/pu_target.GetSumOfWeights());
    pu_source = ref_file.FindObjectAny("hNPUTrue")
    pu_source.Rebin(10)
    pu_source.Scale(1.0/pu_source.GetSumOfWeights())
    pu_weights = pu_target.Clone('pu_weights')
    pu_weights.Divide(pu_source)

    # target file
    new_file_name = output_dir+dataset+'.root'
    if os.path.exists(new_file_name):
        #continue ## FIXME
        if raw_input(new_file_name+' exists, recreate it? [Y/n]')=='n': continue
        
    new_file = TFile(new_file_name,'RECREATE') 
    
    # process the trees
    tree_names = {'MJetTree':'tree','MJetTreeDiLepton':'treeZll','MJetTreeWlnu':'treeWln',}
    for tree_name in tree_names:
        print tree_name
        tree = ref_file.FindObjectAny(tree_name)
        new_tree = tree.CloneTree(-1,'fast')
        new_tree.SetName(tree_names[tree_name])

        # pu weight
        puWeight = array('f',[1.0]) 
        puWeightBranch = new_tree.Branch('puWeight',puWeight,'puWeight/F') 

        # xsec
        xsec = array('f',[1.0]) 
        xsecBranch = new_tree.Branch('xsec',xsec,'xsec/F') 
        xsec[0] = ref_xsec

        # initial number of events
        initialEntries= array('f',[1.0]) 
        initialEntriesBranch = new_tree.Branch('initialEntries',initialEntries,'initialEntries/F')
        initialEntries[0] = ref_initialEntries

        # sf weight
        muonSF = array('f',[1.0]) 
        muonSFBranch = new_tree.Branch('muonSF',muonSF,'muonSF/F') 

        for event in range(0, tree.GetEntries()):
            if event%100000==0: print event 
            tree.GetEntry(event)
            puWeight[0] = pu_weights.GetBinContent(pu_weights.FindFixBin(tree.npu))
            if tree_name != 'tree':
                for bin in muonSFs:
                    if bin[0][0]<=tree.lep1.pt()<=bin[0][1] and bin[1][0]<=abs(tree.lep1.eta())<=bin[1][1]:
                        muonSF[0] = muonSFs[bin]
                if tree_name == 'treeZll':
                    for bin in muonSFs:
                        if bin[0][0]<=tree.lep2.pt()<=bin[0][1] and bin[1][0]<=abs(tree.lep2.eta())<=bin[1][1]:
                            muonSF[0]*=muonSFs[bin]

            puWeightBranch.Fill()
            xsecBranch.Fill()
            initialEntriesBranch.Fill()
            muonSFBranch.Fill()

        new_file.cd()
        new_tree.Write()

    new_file.Close()

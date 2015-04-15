import ROOT
import os
import sys

ROOT.gSystem.Load('libMitPlotsInput.so')
ROOT.gSystem.Load('libMitPlotsPlot.so')

ROOT.gStyle.SetOptStat(0)

# configuration
prodCfg = os.environ['MIT_PROD_CFG']
prodHist = os.environ['MIT_PROD_HIST']
userDir = os.environ['MIT_USER_DIR']
plotsDir = os.environ['MIT_PLOT_DIR'] + '/MonoJet'

# make plots directory
try:
    os.makedirs(plotsDir)
except:
    # do nothing if directory exists
    pass

# initialize input sample list
samples = ROOT.mithep.TaskSamples(prodCfg, prodHist + '/' + prodCfg + '/merged')
samples.ReadFile(userDir + '/config')

# initialize plotting task
plotTask = ROOT.mithep.PlotTask(samples, 19700.)

# set up plot parameters
plotTask.SetHistRanges(0., 800., 0., 0.)
plotTask.SetNBins(40)
plotTask.SetAxisTitles('E_{T}^{miss} (GeV)', 'Number of events')
plotTask.SetPngFileName(plotsDir + '/metStack.pdf')

# run plot task
plotTask.Plot(ROOT.mithep.Stacked, 'MetAnalysisMod/hMet', '', '')

# tweak the canvas appearance    
c1 = ROOT.gROOT.GetListOfCanvases().FindObject('c1')
c1.SetBottomMargin(0.15)
c1.SetTopMargin(0.05)
c1.SetLeftMargin(0.12)
c1.Update()

# make legend (should be part of plot task in the future)
legend = ROOT.TLegend(0.7, 0.7, 0.9, 0.95)
legend.SetFillStyle(0)
legend.SetBorderSize(0)

styles = ROOT.mithep.HistStyles()
styles.PreviousStyle()
garbage = []
for iS in range(samples.NSamples()):
    sample = samples.GetSample(iS)
    if sample.Legend().Data():
        styles.NextStyle()
        dummy = ROOT.TH1F('legend_' + sample.Name().Data(), '', 1, 0., 1.)
        garbage.append(dummy)
        styles.ApplyCurrentStyle(dummy)
        dummy.SetLineWidth(2)
        options = 'FL'
        if dummy.GetFillStyle() == 0:
            options = 'L'
        legend.AddEntry(dummy, sample.Legend().Data(), options)
        
legend.Draw()

c1.Update()

# wait for some input so that we can keep the canvas displayed
sys.stdin.readline()

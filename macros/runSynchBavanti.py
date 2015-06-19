from MitAna.TreeMod.bambu import mithep, analysis

import os

isData = False

goodPVFilterMod = mithep.GoodPVFilterMod(
    MinVertexNTracks = 0,
    MinNDof = 4,
    MaxAbsZ = 24.,
    MaxRho = 2.,
    IsMC = not isData,
    VertexesName = 'PrimaryVertexes'
)

eleIdMod = mithep.ElectronIDMod(
    PtMin = 10.,
    EtaMax = 2.5,
    ApplyEcalFiducial = True,
    IDType = mithep.ElectronTools.kPhys14Veto,
    IsoType = mithep.ElectronTools.kPhys14VetoIso,
    ApplyConversionFilterType1 = True,
    ApplyD0Cut = True,
    ApplyDZCut = True,
    WhichVertex = 0,
    GoodElectronsName = "VetoElectrons",
    RhoAlgo = mithep.PileupEnergyDensity.kFixedGridFastjetAll
)

muonIdIsoMod = mithep.MuonIDMod(
    OutputName = "HWWMuons",
    IntRadius = 0.0,
    ClassType = "GlobalTracker",
    IDType = "WWMuIdV4",
    IsoType = "IsoRingsV0_BDTG_Iso",
    ApplyD0Cut = True,
    D0Cut = 0.02,
    ApplyDZCut = True,
    DZCut = 0.1,
    WhichVertex = 0,
    RhoAlgo = mithep.PileupEnergyDensity.kFixedGridFastjetAll,
    PtMin = 10.,
    EtaCut = 2.4
)

photonIdMod = mithep.PhotonIDMod(
    PtMin = 10.0,
    OutputName = "GoodPhotons",
    IDType = "EgammaMedium",
    IsoType = "NoIso",
    ApplyElectronVeto = True,
    ApplyPixelSeed = False,
    ApplyConversionId = False,
    ApplyFiduciality = True,
    IsData = isData,
    PhotonsFromBranch = True
)

pftauIdMod = mithep.PFTauIDMod(
    PFTausName = "HPSTaus",
    IsLooseId = True,
    IsHPSSel = True,
    PtMin = 15
)
  
#pubAk4Jet = getattr(mithep, 'PublisherMod<mithep::PFJet,mithep::Jet>')('Akt4PFJetsPublisher')
#pubAk4Jet.SetInputName = "AKt4PFJets",
#pubAk4Jet.SetOutputName = "PubAKt4PFJets",
#
#pubAk8Jet = getattr(mithep, 'PublisherMod<mithep::PFJet,mithep::Jet>')('Akt8PFJetsPublisher')
#pubAk8Jet.SetInputName = "AKt8PFJets",
#pubAk8Jet.SetOutputName = "PubAKt8PFJets",
#        
#ak4JetId = mithep.JetIDMod(
#    InputName = pubAk4Jet.GetOutputName(),
#    PtCut = 30.0,
#    EtaMaxCut = 2.5,
#    JetEEMFractionMinCut = 0.00,
#    OutputName = "GoodAk4Jets",
#    ApplyPFLooseId = True
#)
#
#ak8JetId = mithep.JetIDMod(
#    InputName = pubAk8Jet.GetOutputName(),
#    PtCut = 30.0,
#    EtaMaxCut = 2.5,
#    JetEEMFractionMinCut = 0.00,
#    OutputName = "GoodAk8Jets"
#)
#
outMod = mithep.OutputMod(
    UseBrDep = False,
    KeepTamBr = False,
    FileName = 'ntuples.root',
    MaxFileSize = 4096
)
outMod.Drop("*")
outMod.Keep(mithep.Names.gkMCEvtInfoBrn)
outMod.Keep(mithep.Names.gkMCPartBrn)
outMod.Keep(mithep.Names.gkPVBeamSpotBrn)
outMod.Keep(mithep.Names.gkPileupInfoBrn)
outMod.Keep(mithep.Names.gkPileupEnergyDensityBrn)
outMod.Keep("PFMet")

readMC = mithep.ReadMCWeights()

analysis.setSequence(
    goodPVFilterMod * 
    muonIdIsoMod * 
    eleIdMod * 
    photonIdMod * 
    pftauIdMod * 
#    pubAk4Jet * 
#    ak4JetId
    outMod *
    readMC
)

analysis.SetUseMC(True)

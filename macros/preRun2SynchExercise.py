from MitAna.TreeMod.bambu import mithep, analysis
import os

mitdata = os.environ['MIT_DATA']

goodPV = mithep.GoodPVFilterMod(
    MinVertexNTracks = 0,
    MinNDof = 4,
    MaxAbsZ = 24.,
    MaxRho = 2.,
    IsMC = True,
    VertexesName = mithep.Names.gkPVBrn
)

pfPU = mithep.SeparatePileUpMod(
    PFNoPileUpName = "pfNoPU",
    PFPileUpName = "pfPU",
    CheckClosestZVertex = False
)

eleId = mithep.ElectronIdMod(
    OutputName = 'VetoElectrons',
    IdType = mithep.ElectronTools.kPhys14Veto,
    IsoType = mithep.ElectronTools.kPhys14VetoIso,
    ApplyEcalFiducial = True,
    WhichVertex = 0,
    PtMin = 10.,
    EtaMax = 2.5,
    ConversionsName = 'Conversions'
)

muId = mithep.MuonIdMod(
    OutputName = 'LooseMuons',
    IdType = mithep.MuonTools.kNoId,
    IsoType = mithep.MuonTools.kPFIsoBetaPUCorrected,
    PFNoPileupCandidatesName = 'pfNoPU',
    PFPileupCandidatesName = 'pfPU',
    PtMin = 10.,
    EtaMax = 2.4
)

tauId = mithep.PFTauIdMod(
    OutputName = 'LooseTaus',
    PtMin = 18.,
    EtaMax = 2.3
)
tauId.AddDiscriminator(mithep.PFTau.kDiscriminationByDecayModeFindingNewDMs)

photonId = mithep.PhotonIdMod(
    OutputName = 'LoosePhotons',
    IdType = mithep.PhotonTools.kPhys14Loose,
    IsoType = mithep.PhotonTools.kPhys14Loose,
    PtMin = 15.,
    EtaMax = 2.5
)

jetCorr = mithep.JetCorrectionMod(
    InputName = 'AKt4PFJetsCHS',
    CorrectedJetsName = 'AKt4PFJetsCHSL1L2L3',
    RhoAlgo = mithep.PileupEnergyDensity.kFixedGridFastjetAll
)
jetCorr.AddCorrectionFromFile(mitdata + "/MCRUN2_74_V9_L1FastJet_AK4PFchs.txt")
jetCorr.AddCorrectionFromFile(mitdata + "/MCRUN2_74_V9_L2Relative_AK4PFchs.txt")
jetCorr.AddCorrectionFromFile(mitdata + "/MCRUN2_74_V9_L3Absolute_AK4PFchs.txt")

jetId = mithep.JetIdMod(
    InputName = jetCorr.GetOutputName(),
    OutputName = 'CleanedJets',
    UseL1Correction = False,
    UseL2Correction = False,
    UseL3Correction = False,
    ApplyPFLooseId = True,
    ApplyMVACut = True,
    PtMin = 30.,
    EtaMax = 2.5
)

metCorr = mithep.MetCorrectionMod(
    CorrectedName = 'PFType1CorrectedMet',
    JetsName = 'AKt4PFJetsCHS',
    CorrectedJetsName = jetCorr.GetOutputName()
)
metCorr.ApplyType1(True)
metCorr.IsData(False)

preRun2Sync = mithep.PreRun2SynchExercise(
    VerticesName = mithep.Names.gkPVBrn,
    MetName = metCorr.GetOutputName(),
    JetsName = jetId.GetOutputName(),
    ElectronsName = eleId.GetOutputName(),
    MuonsName = muId.GetOutputName(),
    TausName = tauId.GetOutputName(),
    PhotonsName = photonId.GetOutputName()
)

analysis.setSequence(
    goodPV *
    pfPU *
    eleId *
    muId *
    tauId *
    photonId *
    jetCorr *
    jetId *
    metCorr *
    preRun2Sync
)

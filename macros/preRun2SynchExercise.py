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

jetPub = mithep.PFJetToJetPublisherMod(
    Name = 'AKt4PFJetsPublisher',
    InputName = mithep.Names.gkPFJetBrn,
    OutputName = 'AKt4PFJetsDownCasted'
)

jetCorr = mithep.JetCorrectionMod(
    InputName = jetPub.GetOutputName(),
    CorrectedJetsName = 'AKt4PFJetsL1L2L3',
    RhoAlgo = mithep.PileupEnergyDensity.kFixedGridFastjetAll
)
jetCorr.AddCorrectionFromFile(mitdata + "/Summer13_V1_MC_L1FastJet_AK5PF.txt")
jetCorr.AddCorrectionFromFile(mitdata + "/Summer13_V1_MC_L2Relative_AK5PF.txt")
jetCorr.AddCorrectionFromFile(mitdata + "/Summer13_V1_MC_L3Absolute_AK5PF.txt")

jetId = mithep.JetIdMod(
    InputName = jetCorr.GetOutputName(),
    OutputName = 'CleanedJets',
    ApplyPFLooseId = True,
    ApplyMVACut = True,
    PtMin = 30.,
    EtaMax = 2.5
)

metCorr = mithep.MetCorrectionMod(
    CorrectedName = 'PFType1CorrectedMet',
    JetsName = jetPub.GetOutputName(),
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
    eleId *
    muId *
    tauId *
    photonId *
    jetPub *
    jetCorr *
    jetId *
    metCorr *
    preRun2Sync
)

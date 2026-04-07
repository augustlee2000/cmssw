"""
Scouting Hadronic Top DQM Configuration

Self-contained configuration for hadronic top reconstruction from scouting data.
Includes: PF candidate conversion, jet clustering, ParticleNet b-tagging, and DQM analyzer.

Author: August Lee
Date: March 2026
"""

import FWCore.ParameterSet.Config as cms

# ============================================================================
# Step 1: Convert Scouting Particles to Reco PF Candidates
# ============================================================================
scoutingPFCands = cms.EDProducer(
    "Run3ScoutingParticleToRecoPFCandidateProducer",
    scoutingparticle = cms.InputTag("hltScoutingPFPacker"),
)

# ============================================================================
# Step 2: Cluster AK4 Jets from Scouting PF Candidates
# ============================================================================
from RecoJets.JetProducers.ak4PFJets_cfi import ak4PFJets
ak4ScoutingJets = ak4PFJets.clone(
    src = "scoutingPFCands",
    jetPtMin = 20,
)

# ============================================================================
# Step 3: Produce Jet Tag Info for ParticleNet (using scouting track features)
# ============================================================================
ak4ScoutingJetParticleNetJetTagInfos = cms.EDProducer(
    "DeepBoostedJetTagInfoProducer",
    jet_radius = cms.double(0.4),
    min_jet_pt = cms.double(5.0),
    max_jet_eta = cms.double(2.5),
    min_pt_for_track_properties = cms.double(0.95),
    min_pt_for_pfcandidates = cms.double(0.1),
    use_puppiP4 = cms.bool(False),
    include_neutrals = cms.bool(True),
    sort_by_sip2dsig = cms.bool(False),
    min_puppi_wgt = cms.double(-1.0),
    flip_ip_sign = cms.bool(False),
    sip3dSigMax = cms.double(-1.0),
    use_hlt_features = cms.bool(False),
    pf_candidates = cms.InputTag("scoutingPFCands"),
    jets = cms.InputTag("ak4ScoutingJets"),
    puppi_value_map = cms.InputTag(""),
    # Scouting-specific track features from the PF candidate producer
    use_scouting_features = cms.bool(True),
    normchi2_value_map = cms.InputTag("scoutingPFCands", "normchi2"),
    dz_value_map = cms.InputTag("scoutingPFCands", "dz"),
    dxy_value_map = cms.InputTag("scoutingPFCands", "dxy"),
    dzsig_value_map = cms.InputTag("scoutingPFCands", "dzsig"),
    dxysig_value_map = cms.InputTag("scoutingPFCands", "dxysig"),
    lostInnerHits_value_map = cms.InputTag("scoutingPFCands", "lostInnerHits"),
    quality_value_map = cms.InputTag("scoutingPFCands", "quality"),
    trkPt_value_map = cms.InputTag("scoutingPFCands", "trkPt"),
    trkEta_value_map = cms.InputTag("scoutingPFCands", "trkEta"),
    trkPhi_value_map = cms.InputTag("scoutingPFCands", "trkPhi"),
)

# ============================================================================
# Step 4: Run ParticleNet ONNX Model for B-Tagging
# ============================================================================
ak4ScoutingJetParticleNetJetTags = cms.EDProducer(
    "BoostedJetONNXJetTagsProducer",
    jets = cms.InputTag("ak4ScoutingJets"),
    produceValueMap = cms.untracked.bool(True),
    src = cms.InputTag("ak4ScoutingJetParticleNetJetTagInfos"),
    preprocess_json = cms.string("data/preprocess.json"),
    model_path = cms.FileInPath("data/final_model.onnx"),
    flav_names = cms.vstring(["probb", "probc", "probudsg", "probundef"]),
    debugMode = cms.untracked.bool(False),
)

# ============================================================================
# Step 5: Hadronic Top DQM Analyzer
# ============================================================================
scoutingHadronicTopAnalyzer = cms.EDProducer(
    "ScoutingHadronicTopMassAnalyzer",
    OutputPath = cms.string("HLT/ScoutingHadronicTop"),
    
    # Jet inputs (from scouting->reco conversion + jet clustering chain)
    JetCollection = cms.InputTag("ak4ScoutingJets"),
    BTagProbB = cms.InputTag("ak4ScoutingJetParticleNetJetTags", "probb"),
    
    # Muon and MET inputs (native scouting types from file)
    ScoutingMuons = cms.InputTag("hltScoutingMuonPackerVtx"),
    MET = cms.InputTag("hltScoutingPFPacker", "pfMetPt"),
    
    # Trigger configuration
    TriggerResults = cms.InputTag("TriggerResults", "", "HLT"),
    TriggerSelection = cms.vstring(),  # Empty = no trigger filter
    
    # Kinematic selections
    MinJetPt = cms.double(30.0),       # Minimum jet pT [GeV]
    MaxJetEta = cms.double(2.4),       # Maximum jet |eta|
    BTagWP = cms.double(0.5),          # ParticleNet P(b) threshold
    MinMuonPt = cms.double(20.0),      # Minimum muon pT [GeV]
)

# ============================================================================
# Sequences and Tasks
# ============================================================================

# Task containing all producers (for unscheduled execution)
scoutingHadronicTopTask = cms.Task(
    scoutingPFCands,
    ak4ScoutingJets,
    ak4ScoutingJetParticleNetJetTagInfos,
    ak4ScoutingJetParticleNetJetTags,
)

# Sequence for scheduled execution
scoutingHadronicTopSequence = cms.Sequence(
    scoutingPFCands +
    ak4ScoutingJets +
    ak4ScoutingJetParticleNetJetTagInfos +
    ak4ScoutingJetParticleNetJetTags
)

# Full sequence including the analyzer
scoutingHadronicTopDQMSequence = cms.Sequence(
    scoutingHadronicTopSequence +
    scoutingHadronicTopAnalyzer
)


# ============================================================================
# Helper function to customize for different working points
# ============================================================================
def customizeBTagWP(process, wp="medium"):
    """
    Customize b-tagging working point.
    wp: "loose" (0.25), "medium" (0.5), "tight" (0.75)
    """
    wp_values = {"loose": 0.25, "medium": 0.5, "tight": 0.75}
    if wp not in wp_values:
        raise ValueError(f"Unknown working point: {wp}. Use 'loose', 'medium', or 'tight'.")
    
    process.scoutingHadronicTopAnalyzer.BTagWorkingPoint = cms.double(wp_values[wp])
    return process


def customizeTriggerSelection(process, triggers):
    """
    Set trigger selection for the analyzer.
    triggers: list of HLT path substrings, e.g., ["HLT_PFJet500", "HLT_PFHT"]
    """
    process.scoutingHadronicTopAnalyzer.TriggerSelection = cms.vstring(triggers)
    return process

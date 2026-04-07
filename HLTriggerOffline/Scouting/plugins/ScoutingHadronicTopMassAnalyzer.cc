/*
Scouting Hadronic Top Analyzer 
Made by: August Lee University of Notre Dame, alee43@nd.edu
This analyzer reconstructs hadronic and leptonic top quark candidates in scouting data, applying selection cuts and filling histograms for DQM monitoring. It uses jets, muons, MET, and b-tagging information to identify top quark decays.
*/

#include <algorithm>
#include <atomic>
#include <cmath>
#include <string>
#include <vector>

#include "DQMServices/Core/interface/DQMGlobalEDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/Scouting/interface/Run3ScoutingMuon.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "FWCore/Common/interface/TriggerNames.h"


//include statment for the hadronic 2d PDF
#include "data/TopPDF.h"

// Histogram for DQM monitoring
struct SimpleHistos {
  dqm::reco::MonitorElement* hEventCount;
  
  // Jet histograms
  dqm::reco::MonitorElement* hNJets;
  dqm::reco::MonitorElement* hJetPt;
  dqm::reco::MonitorElement* hJetEta;
  dqm::reco::MonitorElement* hJetPhi;
  dqm::reco::MonitorElement* hBTagScore;
  
  // Muon histograms
  dqm::reco::MonitorElement* hNMuons;
  dqm::reco::MonitorElement* hMuonPt;
  dqm::reco::MonitorElement* hMuonEta;
  dqm::reco::MonitorElement* hMuonPhi;
  
  // MET histogram
  dqm::reco::MonitorElement* hMET;

  // Hadronic top histograms (before b-tag)
  dqm::reco::MonitorElement* hHadTopMass;
  dqm::reco::MonitorElement* hHadWMass;

  // Leptonic top histograms (before b-tag)
  dqm::reco::MonitorElement* hLepTopMass;

  // After both b-jets pass b-tag
  dqm::reco::MonitorElement* hHadTopMassBTag;
  dqm::reco::MonitorElement* hHadWMassBTag;
  dqm::reco::MonitorElement* hLepTopMassBTag;
};

class ScoutingHadronicTopMassAnalyzer : public DQMGlobalEDAnalyzer<SimpleHistos> {
public:
  explicit ScoutingHadronicTopMassAnalyzer(const edm::ParameterSet& conf);
  ~ScoutingHadronicTopMassAnalyzer() override;
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void dqmAnalyze(const edm::Event& e, const edm::EventSetup& c, SimpleHistos const&) const override;
  void bookHistograms(DQMStore::IBooker&, edm::Run const&, edm::EventSetup const&, SimpleHistos&) const override;

  // Tokens
  const edm::EDGetTokenT<reco::PFJetCollection> jetsToken_;
  const edm::EDGetTokenT<edm::ValueMap<float>> bTagProbBToken_;
  const edm::EDGetTokenT<std::vector<Run3ScoutingMuon>> muonsToken_;
  const edm::EDGetTokenT<double> metToken_;
  const edm::EDGetTokenT<double> metPhiToken_;
  const edm::EDGetTokenT<edm::TriggerResults> triggerToken_;
  
  // Parameters
  std::string outputPath_;
  std::vector<std::string> triggerSelection_;
  double minJetPt_;
  double maxJetEta_;
  double bTagWP_;
  double minMuonPt_;

  // Counters
  mutable std::atomic<int> nAll_{0};
  mutable std::atomic<int> nTrigger_{0};
  mutable std::atomic<int> nJetsValid_{0};
  mutable std::atomic<int> nBTagValid_{0};
  mutable std::atomic<int> nGe4Jets_{0};
  mutable std::atomic<int> nGe2BTag_{0};
  mutable std::atomic<int> nMuonsValid_{0};
  mutable std::atomic<int> nMuonSel_{0};
  mutable std::atomic<int> nMETCut_{0};
  mutable std::atomic<int> nHadTop_{0};
  mutable std::atomic<int> nLepTop_{0};
  mutable std::atomic<int> nBothBTag_{0};
};

ScoutingHadronicTopMassAnalyzer::ScoutingHadronicTopMassAnalyzer(const edm::ParameterSet& iConfig)
    : jetsToken_(consumes<reco::PFJetCollection>(iConfig.getParameter<edm::InputTag>("JetCollection"))),
      bTagProbBToken_(consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("BTagProbB"))),
      muonsToken_(consumes<std::vector<Run3ScoutingMuon>>(iConfig.getParameter<edm::InputTag>("ScoutingMuons"))),
      metToken_(consumes<double>(iConfig.getParameter<edm::InputTag>("MET"))),
      metPhiToken_(consumes<double>(iConfig.getParameter<edm::InputTag>("METPhi"))),
      triggerToken_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("TriggerResults"))),
      outputPath_(iConfig.getParameter<std::string>("OutputPath")),
      triggerSelection_(iConfig.getParameter<std::vector<std::string>>("TriggerSelection")),
      minJetPt_(iConfig.getParameter<double>("MinJetPt")),
      maxJetEta_(iConfig.getParameter<double>("MaxJetEta")),
      bTagWP_(iConfig.getParameter<double>("BTagWP")),
      minMuonPt_(iConfig.getParameter<double>("MinMuonPt")) {}

ScoutingHadronicTopMassAnalyzer::~ScoutingHadronicTopMassAnalyzer() {
  edm::LogPrint("ScoutingHadronicTopMassAnalyzer")
    << "\n=== CutFlow Summary ==="
    << "\n  All events:         " << nAll_
    << "\n  Pass trigger:       " << nTrigger_
    << "\n  Jets valid:         " << nJetsValid_
    << "\n  BTag valid:         " << nBTagValid_
    << "\n  >= 4 jets:          " << nGe4Jets_
    << "\n  >= 2 b-tags:        " << nGe2BTag_
    << "\n  Muons valid:        " << nMuonsValid_
    << "\n  Muon selection:     " << nMuonSel_
    << "\n  MET cut:            " << nMETCut_
    << "\n  Hadronic top found: " << nHadTop_
    << "\n  Leptonic top found: " << nLepTop_
    << "\n  Both b-tag pass:    " << nBothBTag_
    << "\n========================";
}

void ScoutingHadronicTopMassAnalyzer::bookHistograms(DQMStore::IBooker& ibook, 
                                                      edm::Run const&, 
                                                      edm::EventSetup const&, 
                                                      SimpleHistos& histos) const {
  ibook.setCurrentFolder(outputPath_);
  
  histos.hEventCount = ibook.book1D("EventCount", "Event Count", 1, 0, 1);
  
  // Jet histograms
  histos.hNJets = ibook.book1D("NJets", "Number of Jets", 20, 0, 20);
  histos.hJetPt = ibook.book1D("JetPt", "Jet pT", 100, 0, 500);
  histos.hJetEta = ibook.book1D("JetEta", "Jet Eta", 50, -2.5, 2.5);
  histos.hJetPhi = ibook.book1D("JetPhi", "Jet Phi", 64, -3.14, 3.14);
  histos.hBTagScore = ibook.book1D("BTagScore", "B-Tag Score", 100, 0, 1);
  
  // Muon histograms
  histos.hNMuons = ibook.book1D("NMuons", "Number of Muons", 10, 0, 10);
  histos.hMuonPt = ibook.book1D("MuonPt", "Muon pT", 100, 0, 500);
  histos.hMuonEta = ibook.book1D("MuonEta", "Muon Eta", 50, -2.5, 2.5);
  histos.hMuonPhi = ibook.book1D("MuonPhi", "Muon Phi", 64, -3.14, 3.14);
  
  // MET histogram
  histos.hMET = ibook.book1D("MET", "Missing ET", 100, 0, 500);

  // Hadronic top histograms (before b-tag)
  histos.hHadTopMass = ibook.book1D("HadTopMass", "Hadronic Top Mass (before b-tag)", 100, 0, 500);
  histos.hHadWMass = ibook.book1D("HadWMass", "Hadronic W Mass (before b-tag)", 100, 0, 300);

  // Leptonic top histogram (before b-tag)
  histos.hLepTopMass = ibook.book1D("LepTopMass", "Leptonic Top Mass (before b-tag)", 100, 0, 500);

  // After both b-jets pass b-tag
  histos.hHadTopMassBTag = ibook.book1D("HadTopMassBTag", "Hadronic Top Mass (both b-tag pass)", 100, 0, 500);
  histos.hHadWMassBTag = ibook.book1D("HadWMassBTag", "Hadronic W Mass (both b-tag pass)", 100, 0, 300);
  histos.hLepTopMassBTag = ibook.book1D("LepTopMassBTag", "Leptonic Top Mass (both b-tag pass)", 100, 0, 500);
}

void ScoutingHadronicTopMassAnalyzer::dqmAnalyze(const edm::Event& iEvent,
                                                   const edm::EventSetup& iSetup,
                                                   SimpleHistos const& histos) const {
  histos.hEventCount->Fill(0.5);
  ++nAll_;
  
  // Check trigger if configured
  if (!triggerSelection_.empty()) {
    edm::Handle<edm::TriggerResults> triggerResults;
    iEvent.getByToken(triggerToken_, triggerResults);
    
    if (!triggerResults.isValid())
      return;
    
    const edm::TriggerNames& names = iEvent.triggerNames(*triggerResults);
    bool passTrigger = false;
    
    for (unsigned int i = 0; i < triggerResults->size(); ++i) {
      if (!triggerResults->accept(i)) continue;
      
      const std::string& name = names.triggerName(i);
      for (const auto& sel : triggerSelection_) {
        if (name.find(sel) != std::string::npos) {
          passTrigger = true;
          break;
        }
      }
      if (passTrigger) break;
    }
    
    if (!passTrigger) return;
  }
  ++nTrigger_;
  
  // Get reco jets (produced by scouting→reco conversion chain)
  edm::Handle<reco::PFJetCollection> jets;
  iEvent.getByToken(jetsToken_, jets);
  if (!jets.isValid()) {
    edm::LogWarning("ScoutingHadronicTopMassAnalyzer") 
      << "Jets collection not found in event " << iEvent.id().event();
    return;
  }
  ++nJetsValid_;
  
  // Get b-tag scores (ParticleNet probb ValueMap keyed to ak4ScoutingJets)
  edm::Handle<edm::ValueMap<float>> bTagProbB;
  iEvent.getByToken(bTagProbBToken_, bTagProbB);
  bool hasBTag = bTagProbB.isValid();
  if (!hasBTag) {
    edm::LogWarning("ScoutingHadronicTopMassAnalyzer") 
      << "B-tag scores not found in event " << iEvent.id().event();
  } else {
    ++nBTagValid_;
  }
  
  // Count and fill jet histograms
  int nJets = 0;
  int nBTagJets = 0;
  edm::LogPrint("ScoutingHadronicTopMassAnalyzer") 
    << "Event " << iEvent.id().event() << " - B-Tag Scores (hasBTag=" << hasBTag << "):";
  
  for (size_t i = 0; i < jets->size(); ++i) {
    const auto& jet = (*jets)[i];
    
    if (jet.pt() > minJetPt_ && std::abs(jet.eta()) < maxJetEta_) {
      nJets++;
      histos.hJetPt->Fill(jet.pt());
      histos.hJetEta->Fill(jet.eta());
      histos.hJetPhi->Fill(jet.phi());
      
      if (hasBTag) {
        edm::Ref<reco::PFJetCollection> jetRef(jets, i);
        float probB = (*bTagProbB)[jetRef];
        histos.hBTagScore->Fill(probB);
        if (probB > bTagWP_)
          nBTagJets++;
        edm::LogPrint("ScoutingHadronicTopMassAnalyzer")
          << "  Jet " << (nJets-1) << ": pT=" << jet.pt()
          << " eta=" << jet.eta()
          << " probb=" << probB;
      }
    }
  }
  histos.hNJets->Fill(nJets);
  
  // Require at least 4 jets passing cuts
  if (nJets < 4) {
    return;
  }
  ++nGe4Jets_;
  
  // Require at least 2 b-tagged jets
  if (nBTagJets < 2) {
    return;
  }
  ++nGe2BTag_;
  
  // Get scouting muons
  edm::Handle<std::vector<Run3ScoutingMuon>> muons;
  iEvent.getByToken(muonsToken_, muons);
  if (!muons.isValid()) {
    edm::LogWarning("ScoutingHadronicTopMassAnalyzer") 
      << "Muons collection not found in event " << iEvent.id().event();
    return;
  }
  ++nMuonsValid_;
  
  // Require exactly one muon above 30 GeV and all others below 15 GeV
  int nMuonsAbove30 = 0;
  int nMuonsBelow15 = 0;
  int nMuons = 0;
  
  for (const auto& muon : *muons) {
    if (muon.pt() > 30.0) {
      nMuonsAbove30++;
    } else if (muon.pt() < 15.0) {
      nMuonsBelow15++;
    } else {
      // Muon is between 15 and 30 GeV, does not pass selection
      return;
    }
  }
  
  // Check the constraint: exactly one above 30 GeV, rest below 15 GeV
  if (nMuonsAbove30 != 1) {
    return;
  }
  ++nMuonSel_;
  
  nMuons = nMuonsAbove30 + nMuonsBelow15;
  
  // Fill muon histograms for all muons passing the selection
  for (const auto& muon : *muons) {
    if (muon.pt() > 30.0 || muon.pt() < 15.0) {
      histos.hMuonPt->Fill(muon.pt());
      histos.hMuonEta->Fill(muon.eta());
      histos.hMuonPhi->Fill(muon.phi());
    }
  }
  histos.hNMuons->Fill(nMuons);
  
  // Find the high-pT muon (above 30 GeV) for MET cut
  double highPtMuonPt = -1.0;
  for (const auto& muon : *muons) {
    if (muon.pt() > 30.0) {
      highPtMuonPt = muon.pt();
      break;
    }
  }
  
  // Get MET and apply linear cut: MET > -2 * muon_pt + 80
  edm::Handle<double> met;
  iEvent.getByToken(metToken_, met);
  if (!met.isValid())
    return;
  double metCutThreshold = -2.0 * highPtMuonPt + 80.0;
  if (*met <= metCutThreshold)
    return;
  ++nMETCut_;
  histos.hMET->Fill(*met);

  // --- Find best hadronic top: loop all 3-jet combos, pick highest 2D PDF score ---
  // Collect the indices (into the full jets collection) of jets passing cuts
  std::vector<size_t> goodJetIdx;
  for (size_t i = 0; i < jets->size(); ++i) {
    const auto& jet = (*jets)[i];
    if (jet.pt() > minJetPt_ && std::abs(jet.eta()) < maxJetEta_)
      goodJetIdx.push_back(i);
  }

  double bestProb = -1.0;
  double bestTopMass = -1.0;
  double bestWMass = -1.0;
  int bestHadB = -1, bestW1 = -1, bestW2 = -1;

  int nGood = static_cast<int>(goodJetIdx.size());
  for (int ib = 0; ib < nGood; ++ib) {
    for (int iw1 = 0; iw1 < nGood; ++iw1) {
      if (iw1 == ib) continue;
      for (int iw2 = iw1 + 1; iw2 < nGood; ++iw2) {
        if (iw2 == ib) continue;

        const auto& jb  = (*jets)[goodJetIdx[ib]];
        const auto& jw1 = (*jets)[goodJetIdx[iw1]];
        const auto& jw2 = (*jets)[goodJetIdx[iw2]];

        auto pW  = jw1.p4() + jw2.p4();
        auto pTop = pW + jb.p4();
        double wMass   = pW.mass();
        double topMass = pTop.mass();

        double prob = TopPDFs::hadronicPdf(topMass, wMass);
        if (prob > bestProb) {
          bestProb    = prob;
          bestTopMass = topMass;
          bestWMass   = wMass;
          bestHadB    = ib;
          bestW1      = iw1;
          bestW2      = iw2;
        }
      }
    }
  }

  if (bestProb > 0) {
    ++nHadTop_;
    histos.hHadTopMass->Fill(bestTopMass);
    histos.hHadWMass->Fill(bestWMass);
    edm::LogPrint("ScoutingHadronicTopMassAnalyzer")
      << "  Best combo: hadB=" << bestHadB << " w1=" << bestW1 << " w2=" << bestW2
      << " topM=" << bestTopMass << " wM=" << bestWMass << " prob=" << bestProb;

    // --- Leptonic top: highest-pT remaining jet + muon + neutrino ---
    // Find highest-pT good jet not used in the hadronic combo
    int bestLepB = -1;
    double bestLepBPt = -1.0;
    for (int ij = 0; ij < nGood; ++ij) {
      if (ij == bestHadB || ij == bestW1 || ij == bestW2) continue;
      double pt = (*jets)[goodJetIdx[ij]].pt();
      if (pt > bestLepBPt) {
        bestLepBPt = pt;
        bestLepB = ij;
      }
    }

    if (bestLepB >= 0) {
      const auto& lepBJet = (*jets)[goodJetIdx[bestLepB]];

      // Get the selected muon (the one above 30 GeV)
      math::PtEtaPhiMLorentzVector muP4;
      for (const auto& muon : *muons) {
        if (muon.pt() > 30.0) {
          muP4 = math::PtEtaPhiMLorentzVector(muon.pt(), muon.eta(), muon.phi(), 0.1057);
          break;
        }
      }

      // Get MET phi for neutrino
      edm::Handle<double> metPhi;
      iEvent.getByToken(metPhiToken_, metPhi);
      if (metPhi.isValid()) {
        double metPt = *met;
        double mPhi  = *metPhi;
        double nu_px = metPt * std::cos(mPhi);
        double nu_py = metPt * std::sin(mPhi);

        // Solve for neutrino pz using W mass constraint
        constexpr double mW = 80.4;
        double mu_px = muP4.px();
        double mu_py = muP4.py();
        double mu_pz = muP4.pz();
        double mu_E  = muP4.energy();

        double A = (mW * mW) / 2.0 + mu_px * nu_px + mu_py * nu_py;
        double mu_pt2 = mu_px * mu_px + mu_py * mu_py;
        double disc = A * A * mu_pz * mu_pz - mu_pt2 * (mu_E * mu_E * metPt * metPt - A * A);

        double nu_pz = 0.0;
        if (disc >= 0) {
          double sol1 = (A * mu_pz + std::sqrt(disc)) / mu_pt2;
          double sol2 = (A * mu_pz - std::sqrt(disc)) / mu_pt2;
          nu_pz = (std::abs(sol1) < std::abs(sol2)) ? sol1 : sol2;
        } else {
          // Complex: take real part
          nu_pz = A * mu_pz / mu_pt2;
        }

        double nu_E = std::sqrt(metPt * metPt + nu_pz * nu_pz);
        math::XYZTLorentzVector nuP4(nu_px, nu_py, nu_pz, nu_E);

        auto lepTopP4 = lepBJet.p4() + math::XYZTLorentzVector(muP4.px(), muP4.py(), muP4.pz(), muP4.energy()) + nuP4;
        double lepTopMass = lepTopP4.mass();
        ++nLepTop_;

        // Fill leptonic top mass before b-tag cut
        histos.hLepTopMass->Fill(lepTopMass);

        edm::LogPrint("ScoutingHadronicTopMassAnalyzer")
          << "  LepTop: lepB=" << bestLepB << " nu_pz=" << nu_pz
          << " lepTopM=" << lepTopMass;

        // Check if both hadronic and leptonic b-jets pass b-tag > 0.208
        if (hasBTag) {
          edm::Ref<reco::PFJetCollection> hadBRef(jets, goodJetIdx[bestHadB]);
          edm::Ref<reco::PFJetCollection> lepBRef(jets, goodJetIdx[bestLepB]);
          float hadBScore = (*bTagProbB)[hadBRef];
          float lepBScore = (*bTagProbB)[lepBRef];
          edm::LogPrint("ScoutingHadronicTopMassAnalyzer")
            << "  HadB b-tag score=" << hadBScore << " LepB b-tag score=" << lepBScore;
          if (hadBScore > 0.208 && lepBScore > 0.208) {
            ++nBothBTag_;
            histos.hHadTopMassBTag->Fill(bestTopMass);
            histos.hHadWMassBTag->Fill(bestWMass);
            histos.hLepTopMassBTag->Fill(lepTopMass);
          }
        }
      }
    }
  }
}

void ScoutingHadronicTopMassAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<std::string>("OutputPath", "HLT/ScoutingHadronicTop");
  desc.add<edm::InputTag>("JetCollection", edm::InputTag("ak4ScoutingJets"));
  desc.add<edm::InputTag>("BTagProbB", edm::InputTag("ak4ScoutingJetParticleNetJetTags", "probb"));
  desc.add<edm::InputTag>("ScoutingMuons", edm::InputTag("hltScoutingMuonPackerVtx"));
  desc.add<edm::InputTag>("MET", edm::InputTag("hltScoutingPFPacker", "pfMetPt"));
  desc.add<edm::InputTag>("METPhi", edm::InputTag("hltScoutingPFPacker", "pfMetPhi"));
  desc.add<edm::InputTag>("TriggerResults", edm::InputTag("TriggerResults", "", "HLT"));
  desc.add<std::vector<std::string>>("TriggerSelection", {});
  desc.add<double>("MinJetPt", 30.0);
  desc.add<double>("MaxJetEta", 2.4);
  desc.add<double>("BTagWP", 0.208);
  desc.add<double>("MinMuonPt", 20.0);
  descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(ScoutingHadronicTopMassAnalyzer);

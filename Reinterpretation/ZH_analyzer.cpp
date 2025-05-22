#include "SampleAnalyzer/User/Analyzer/user.h"
using namespace MA5;
using namespace std;
#include <cstdlib>
#include <Eigen/Dense>
#include <tuple>

const std::vector<std::string> regionsForHistos = {"SR", "SR_TRACKUP", "CRDY", "CRDY_TRACKUP", "CRTT", "CRTT_TRACKUP"};
const std::vector<std::string> regionsForHistos_CENTRAL = {"SR", "CRDY", "CRTT"};
const std::vector<std::string> regionsForHistos_TRACKUP = {"SR_TRACKUP", "CRDY_TRACKUP", "CRTT_TRACKUP"};

void user::FillHistoAllRegions(const std::string &baseName, double value, const std::vector<std::string> &regions)
{
  for (const auto &r : regions) {
    Manager()->FillHisto(baseName + "_" + r, value);
  }
}

template <typename T>
T normalizePhi(T phi)
{
  while (phi > T(M_PI))
    phi -= T(2 * M_PI);
  while (phi <= -T(M_PI))
    phi += T(2 * M_PI);
  return phi;
}

template <typename T, typename U>
T computeDeltaPhi(T phi1, U phi2)
{
  phi1 = normalizePhi(phi1);
  phi2 = normalizePhi(phi2);
  T dphi = phi1 - phi2;
  while (dphi > T(M_PI))
    dphi -= T(2 * M_PI);
  while (dphi <= -T(M_PI))
    dphi += T(2 * M_PI);
  return dphi;
}

template <typename T, typename U>
double calculateDeltaR(const T &obj1, const U &obj2)
{
  double eta1 = obj1.eta();
  double eta2 = obj2.eta();

  double dphi = std::fabs(computeDeltaPhi(obj1.phi(), obj2.phi()));
  double deta = eta1 - eta2;

  return std::sqrt(deta * deta + dphi * dphi);
}

// Returns the mass (in GeV) for a charged particle based on its PDG id
// Common charged leptons, mesons, and baryons are included
// If the PDG id is not recognized, defaults to the charged pion mass
double getMassFromPDG(int pdgid)
{
  int absid = std::abs(pdgid);
  switch (absid)
  {
  case 11:
    return 0.0005109989461; // electron
  case 13:
    return 0.1056583745; // muon
  case 211:
    return 0.13957039; // charged pion
  case 321:
    return 0.493677; // charged kaon
  case 213:
    return 0.77526; // charged ρ meson
  case 323:
    return 0.89166; // K*(892)+ meson
  case 2212:
    return 0.9382720813; // proton
  case 2214:
    return 1.232; // Δ+ baryon
  case 2224:
    return 1.232; // Δ++ baryon
  case 411:
    return 1.86965; // D+ meson
  case 431:
    return 1.96834; // D_s+ meson
  case 3222:
    return 1.18937; // Σ+ baryon
  case 3112:
    return 1.19745; // Σ- baryon
  case 3312:
    return 1.32171; // Ξ- baryon
  case 3334:
    return 1.67245; // Ω- baryon
  case 521:
    return 5.279; // B+ meson
  case 4122:
    return 2.28646; // Λ_c+ baryon
  case 4222:
    return 2.453; // Σ_c++ baryon
  default:
  {
    std::cout << "[FATAL] Unrecognized PDG id in getMassFromPDG: " + std::to_string(pdgid) + " - please update the switch function in getMassFromPDG with this particle's PDG ID and mass (Gev)." << std::endl;
    throw std::runtime_error("[FATAL] Unrecognized PDG id in getMassFromPDG: " + std::to_string(pdgid) + " - please update the switch function in getMassFromPDG with this particle's PDG ID and mass (Gev).");
  }
  }
}

std::vector<fastjet::PseudoJet> boostToSUEP(const std::vector<fastjet::PseudoJet> &constituents, const fastjet::PseudoJet &jet)
{
  std::vector<fastjet::PseudoJet> boosted_particles;

  // Jet four momenta components
  double E_jet = jet.E();
  double px_jet = jet.px();
  double py_jet = jet.py();
  double pz_jet = jet.pz();

  if (E_jet == 0.0)
  {
    throw std::runtime_error("[FATAL] Division by zero in boostToSUEP: E_jet is zero.");
  }

  // Jet boosts
  double bx = px_jet / E_jet;
  double by = py_jet / E_jet;
  double bz = pz_jet / E_jet;

  double beta2 = bx * bx + by * by + bz * bz;
  double gamma = 1.0 / std::sqrt(1.0 - beta2);

  if (beta2 >= 1.0)
  {
    throw std::runtime_error("[FATAL] Invalid boost in boostToSUEP: beta^2 >= 1 would cause sqrt of negative number.");
  }

  // Boosting
  for (const auto &p : constituents)
  {
    double E = p.E();
    double px = p.px();
    double py = p.py();
    double pz = p.pz();

    double pb = px * bx + py * by + pz * bz;

    double E_prime = gamma * (E - pb);
    double coeff = ((gamma - 1.0) * pb / beta2) - gamma * E;

    double px_prime = px + coeff * bx;
    double py_prime = py + coeff * by;
    double pz_prime = pz + coeff * bz;

    fastjet::PseudoJet boosted(px_prime, py_prime, pz_prime, E_prime);
    boosted_particles.push_back(boosted);
  }

  return boosted_particles;
}

double sphericity(const std::vector<fastjet::PseudoJet> &particles, double r)
{

  // Initialize sphericity matrix
  double S_xx = 0.0, S_xy = 0.0, S_xz = 0.0;
  double S_yy = 0.0, S_yz = 0.0, S_zz = 0.0;
  double norm = 0.0;

  // Calculate momentum components & and normalization factor
  for (const auto &particle : particles)
  {
    double px = particle.px();
    double py = particle.py();
    double pz = particle.pz();
    double p = std::sqrt(px * px + py * py + pz * pz);

    double weight = std::pow(p, r - 2.0); // Weight for each component

    // Calculating the sphericitiy tensor components
    S_xx += px * px * weight;
    S_xy += px * py * weight;
    S_xz += px * pz * weight;
    S_yy += py * py * weight;
    S_yz += py * pz * weight;
    S_zz += pz * pz * weight;

    // Normalization
    norm += std::pow(p, r);
  }

  // Normalize the matrix if norm > 0
  if (norm > 0)
  {
    S_xx /= norm;
    S_xy /= norm;
    S_xz /= norm;
    S_yy /= norm;
    S_yz /= norm;
    S_zz /= norm;
  }

  // Construct sphericity matrix
  Eigen::Matrix3d S;
  S << S_xx, S_xy, S_xz,
      S_xy, S_yy, S_yz,
      S_xz, S_yz, S_zz;

  // Get sphericity matrix eigenvalues
  Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> solver(S);
  if (solver.info() != Eigen::Success)
  {
    throw std::runtime_error("[FATAL] Sphericity eigenvalue solver failed — invalid or non-symmetric matrix input.");
  }

  // Sort the eigenvalues
  Eigen::Vector3d eigenvals = solver.eigenvalues();
  std::sort(eigenvals.data(), eigenvals.data() + eigenvals.size());

  // Grab the two smallest eigenvalues
  double eigenval1 = eigenvals[0];
  double eigenval2 = eigenvals[1];

  // Calculate the sphericity
  double sphericityval = 1.5 * (eigenval1 + eigenval2);

  return sphericityval;
}

bool compareBypT(const RecLeptonFormat &a, const RecLeptonFormat &b)
{
  return a.pt() > b.pt(); // Sort by pT in descending order
}

// Function to sort various lepton collections by pT
void sortLeptonCollections(std::vector<RecLeptonFormat> &electrons,
                           std::vector<RecLeptonFormat> &posElectrons,
                           std::vector<RecLeptonFormat> &negElectrons,
                           std::vector<RecLeptonFormat> &muons,
                           std::vector<RecLeptonFormat> &posMuons,
                           std::vector<RecLeptonFormat> &negMuons,
                           std::vector<RecLeptonFormat> &leptons,
                           std::vector<RecLeptonFormat> &posLeptons,
                           std::vector<RecLeptonFormat> &negLeptons)
{
  // Sorting electrons
  std::sort(electrons.begin(), electrons.end(), compareBypT);
  std::sort(posElectrons.begin(), posElectrons.end(), compareBypT);
  std::sort(negElectrons.begin(), negElectrons.end(), compareBypT);

  // Sorting muons
  std::sort(muons.begin(), muons.end(), compareBypT);
  std::sort(posMuons.begin(), posMuons.end(), compareBypT);
  std::sort(negMuons.begin(), negMuons.end(), compareBypT);

  // Sorting leptons
  std::sort(leptons.begin(), leptons.end(), compareBypT);
  std::sort(posLeptons.begin(), posLeptons.end(), compareBypT);
  std::sort(negLeptons.begin(), negLeptons.end(), compareBypT);
}

bool iso(const RecLeptonFormat &lepton,
         const std::vector<RecTrackFormat> &centralTracks,
         const std::vector<RecParticleFormat> &eflowPhotons,
         const std::vector<RecParticleFormat> &eflowNeutralHadrons,
         double iso_minpt,
         double deltaRmax,
         double deltaRmin,
         std::string iso = "")
{

  double ratiomax = 0.0;
  double abs_eta = std::fabs(lepton.eta());

  if (iso == "WP90")
  {
    if (abs_eta < 0.8)
      ratiomax = 0.23;
    else if (abs_eta < 1.44)
      ratiomax = 0.25;
    else // abs_eta >= 1.57
      ratiomax = 0.25;
  }
  else if (iso == "WP80")
  {
    if (abs_eta < 0.8)
      ratiomax = 0.09;
    else if (abs_eta < 1.44)
      ratiomax = 0.09;
    else // abs_eta >= 1.57
      ratiomax = 0.105;
  }
  else if (iso == "pfIso2")
  {
    ratiomax = 0.25;
  }
  else if (iso == "pfIso5")
  {
    ratiomax = 0.1;
  }

  double lep_pt = lepton.pt();
  double totalpt = 0.0;
  // Loop over reconstructed tracks
  for (const auto &track : centralTracks)
  {
    double dr = lepton.dr(track);
    double pt = track.pt();
    double pdgid = fabs(track.pdgid());

    if (pdgid == 11 or pdgid == 13)
    {
      continue;
    }
    if (dr < deltaRmax && dr > deltaRmin && pt > iso_minpt)
      totalpt += pt;
  }

  for (const auto &photon : eflowPhotons)
  {
    double dr = lepton.dr(photon);
    double pt = photon.pt();
    if (dr < deltaRmax && dr > deltaRmin && pt > iso_minpt)
    {
      totalpt += pt;
    }
  }

  for (const auto &neutral : eflowNeutralHadrons)
  {
    double dr = lepton.dr(neutral);
    double pt = neutral.pt();
    if (dr < 0.01)
    {
      continue;
    }
    if (dr < deltaRmax && dr > deltaRmin && pt > iso_minpt)
    {
      totalpt += pt;
    }
  }

  double pt_ratio = totalpt / lep_pt;

  return (pt_ratio <= ratiomax);
}

vector<RecLeptonFormat> filter_muons(vector<RecLeptonFormat> objects,
                                     float ptmin,
                                     float etamax,
                                     float d0,
                                     float dz,
                                     float iso_pTMin,
                                     float iso_dRMax,
                                     float iso_dRMin,
                                     const vector<RecTrackFormat> &centralTracks,
                                     const vector<RecParticleFormat> &eflowPhotons,
                                     const vector<RecParticleFormat> &eflowNeutralHadrons,
                                     std::string charge = "")
{

  // Helper function to select muons
  vector<RecLeptonFormat> filtered;
  for (auto &obj : objects)
  {
    // Charge filter
    if (charge == "+")
    {
      if (obj.charge() > 0)
      {
        filtered.push_back(obj);
      }
      continue;
    }
    else if (charge == "-")
    {
      if (obj.charge() < 0)
      {
        filtered.push_back(obj);
      }
      continue;
    }

    // Object Selections
    if (obj.pt() < ptmin)
      continue;

    if (fabs(obj.eta()) > etamax)
      continue;

    if (fabs(obj.d0()) > d0)
      continue;

    if (fabs(obj.dz()) > dz)
      continue;

    if (!iso(obj, centralTracks, eflowPhotons, eflowNeutralHadrons, iso_pTMin, iso_dRMax, iso_dRMin, "pfIso2"))
    {
      continue;
    }

    filtered.push_back(obj);
  }

  return filtered;
}

vector<RecLeptonFormat> filter_electrons(vector<RecLeptonFormat> objects,
                                         float ptmin,
                                         float etamax,
                                         float iso_pTMin,
                                         float iso_dRMax,
                                         float iso_dRMin,
                                         const vector<RecTrackFormat> &centralTracks,
                                         const vector<RecParticleFormat> &eflowPhotons,
                                         const vector<RecParticleFormat> &eflowNeutralHadrons,
                                         std::string charge = "")
{
  // Helper function to select electrons
  vector<RecLeptonFormat> filtered;

  for (auto &obj : objects)
  {
    // Charge filter
    if (charge == "+")
    {
      if (obj.charge() > 0)
      {
        filtered.push_back(obj);
      }
      continue;
    }
    else if (charge == "-")
    {
      if (obj.charge() < 0)
      {
        filtered.push_back(obj);
      }
      continue;
    }

    // Object Selections
    if (obj.pt() < ptmin)
      continue;

    if (fabs(obj.eta()) > etamax)
      continue;

    // Exclude crack region: electrons with 1.444 < |eta| < 1.566
    if (fabs(obj.eta()) > 1.444 && fabs(obj.eta()) < 1.566)
      continue;

    // Apply impact parameter cuts based on eta
    if (fabs(obj.d0()) > (0.05 + 0.05 * (fabs(obj.eta()) > 1.479)))
      continue;

    if (fabs(obj.dz()) > (0.1 + 0.1 * (fabs(obj.eta()) > 1.479)))
      continue;

    if (!iso(obj, centralTracks, eflowPhotons, eflowNeutralHadrons, iso_pTMin, iso_dRMax, iso_dRMin, "WP90"))
    {
      continue;
    }

    filtered.push_back(obj);
  }

  return filtered;
}

std::pair<std::vector<fastjet::PseudoJet>, std::vector<std::vector<fastjet::PseudoJet>>>
getAk4jets(const std::vector<RecTrackFormat> &EFlowTracks,
           const std::vector<RecParticleFormat> &EFlowPhotons,
           const std::vector<RecParticleFormat> &EFlowNeutralHadrons,
           const std::vector<RecLeptonFormat> &leptons,
           float leptonCleaningDeltaR)
{
  std::vector<fastjet::PseudoJet> input_particles;

  auto passesCleaning = [leptonCleaningDeltaR, &leptons](const fastjet::PseudoJet &particle) -> bool
  {
    for (const auto &lepton : leptons)
    {
      if (calculateDeltaR(lepton, particle) < leptonCleaningDeltaR)
        return false;
    }
    return true;
  };

  for (const auto &track : EFlowTracks)
  {
    fastjet::PseudoJet particle(track.px(), track.py(), track.pz(), track.e());
    if (passesCleaning(particle))
      input_particles.emplace_back(particle);
  }

  for (const auto &photon : EFlowPhotons)
  {
    fastjet::PseudoJet particle(photon.px(), photon.py(), photon.pz(), photon.e());
    if (passesCleaning(particle))
      input_particles.emplace_back(particle);
  }

  for (const auto &hadron : EFlowNeutralHadrons)
  {
    fastjet::PseudoJet particle(hadron.px(), hadron.py(), hadron.pz(), hadron.e());
    if (passesCleaning(particle))
      input_particles.emplace_back(particle);
  }

  fastjet::JetDefinition jet_def(fastjet::antikt_algorithm, 0.4);
  fastjet::ClusterSequence cluster_seq(input_particles, jet_def);

  std::vector<fastjet::PseudoJet> inclusive_jets = sorted_by_pt(cluster_seq.inclusive_jets(15.0));

  std::vector<std::vector<fastjet::PseudoJet>> cluster_constituents;
  for (const auto &jet : inclusive_jets)
  {
    cluster_constituents.push_back(jet.constituents());
  }

  return std::make_pair(inclusive_jets, cluster_constituents);
}

std::pair<std::vector<fastjet::PseudoJet>, std::vector<std::vector<fastjet::PseudoJet>>>
filter_Ak4jetsAndConstituents(const std::vector<fastjet::PseudoJet> &jets,
                              const std::vector<std::vector<fastjet::PseudoJet>> &constituents,
                              float ptmin,
                              float etamax,
                              float deltaRmin,
                              const std::vector<RecLeptonFormat> &leptons)
{
  std::vector<fastjet::PseudoJet> filteredJets;
  std::vector<std::vector<fastjet::PseudoJet>> filteredConstituents;

  for (std::size_t i = 0; i < jets.size(); ++i)
  {
    const auto &jet = jets[i];

    if (jet.pt() < ptmin)
      continue;
    if (std::fabs(jet.eta()) > etamax)
      continue;

    bool passesDeltaRCut = true;
    for (const auto &lepton : leptons)
    {
      if (calculateDeltaR(lepton, jet) <= deltaRmin)
      {
        passesDeltaRCut = false;
        break;
      }
    }

    if (!passesDeltaRCut)
      continue;

    filteredJets.push_back(jet);
    filteredConstituents.push_back(constituents[i]);
  }

  return std::make_pair(filteredJets, filteredConstituents);
}

std::vector<fastjet::PseudoJet> applyEnergyScale(const std::vector<fastjet::PseudoJet> &jets)
{
  std::vector<fastjet::PseudoJet> scaledJets;
  for (auto jet : jets)
  {
    double pt = jet.pt();
    double eta = jet.eta();

    // Compute the scale factor using the given formula:
    // scale = sqrt( (2.5 - 0.15 * |eta|)^2 / pt + 1.5 )
    double scale = sqrt(pow(2.5 - 0.15 * fabs(eta), 2.0) / pt + 1.5);

    // Only apply the scale if it's positive
    if (scale > 0.0)
    {
      jet.reset_momentum(jet.px() * scale, jet.py() * scale, jet.pz() * scale, jet.e() * scale);
    }
    scaledJets.push_back(jet);
  }
  return scaledJets;
}

std::vector<double> computeSoftLeptonMistagSF(const std::vector<fastjet::PseudoJet> &jets, const std::vector<RecTrackFormat> &tracks, double drMatch) {
  std::vector<double> sfVec;
  sfVec.reserve(jets.size());
  for (const auto &jet : jets) {
    double sf = 1.0;  // default: no qualifying soft lepton
    for (const auto &trk : tracks) {
      int absId = std::abs(trk.pdgid());
      if (absId != 11 && absId != 13) continue;
      if (calculateDeltaR(trk, jet) < drMatch)
      {
        if ((absId == 13 || absId == 11) && trk.pt() > 5.0)
        {
          sf = 2.0;
          break;                            
        }
      }
    }
    sfVec.push_back(sf);
  }
  return sfVec;
}

bool BTagVeto(const std::vector<fastjet::PseudoJet> &jets, const std::vector<double> &mistagSF)
{
  // Iterate over jets
  for (std::size_t i = 0; i < jets.size(); ++i)
  {
    const auto &jet = jets[i];
    double pt = jet.pt();
    double abs_eta = fabs(jet.eta());
    double eff = 0.0;

    // Determine the mis-id rate (b-tag efficiency) using the provided piecewise function
    if (pt <= 40.0)
    {
      if (abs_eta <= 1.5)
        eff = 0.144;
      else if (abs_eta <= 2.5)
        eff = 0.236;
    }
    else if (pt <= 50.0)
    {
      if (abs_eta <= 1.5)
        eff = 0.104;
      else if (abs_eta <= 2.5)
        eff = 0.164;
    }
    else if (pt <= 60.0)
    {
      if (abs_eta <= 1.5)
        eff = 0.086;
      else if (abs_eta <= 2.5)
        eff = 0.132;
    }
    else if (pt <= 80.0)
    {
      if (abs_eta <= 1.5)
        eff = 0.075;
      else if (abs_eta <= 2.5)
        eff = 0.114;
    }
    else if (pt <= 100.0)
    {
      if (abs_eta <= 1.5)
        eff = 0.067;
      else if (abs_eta <= 2.5)
        eff = 0.098;
    }
    else if (pt <= 150.0)
    {
      if (abs_eta <= 1.5)
        eff = 0.063;
      else if (abs_eta <= 2.5)
        eff = 0.094;
    }
    else if (pt <= 200.0)
    {
      if (abs_eta <= 1.5)
        eff = 0.068;
      else if (abs_eta <= 2.5)
        eff = 0.108;
    }
    else if (pt <= 500.0)
    {
      if (abs_eta <= 1.5)
        eff = 0.105;
      else if (abs_eta <= 2.5)
        eff = 0.176;
    }
    else // pt > 500.0
    {
      if (abs_eta <= 1.5)
        eff = 0.282;
      else if (abs_eta <= 2.5)
        eff = 0.421;
    }

    // --- apply per‑jet soft‑lepton mistag scale factor ---
    if (i < mistagSF.size())
    {
      eff = std::min(1.0, eff * mistagSF[i]); // cap at 100 %
    }

    // Simulate b-tag decision via a random number.
    // If the random number is less than the efficiency, treat jet as b-tagged.
    double rnd = static_cast<double>(rand()) / RAND_MAX;
    if (rnd < eff)
    {
      return false;
    }
  }
  return true;
}

std::tuple<std::vector<RecTrackFormat>, std::vector<RecTrackFormat>> doTracksDropping(const std::vector<RecTrackFormat> &tracks) {
  const double probLowPt  = 0.03; // low-pt drop probability
  const double probHighPt = 0.01;  // high-pt drop probability

  // Central is just the original tracks
  std::vector<RecTrackFormat> central = tracks;
  std::vector<RecTrackFormat> up;

  // Apply dropping per-track
  for (const auto &trk : tracks) {
    double pt = trk.pt();
    double prob = (pt < 20.0 ? probLowPt : probHighPt);
    double r = static_cast<double>(rand()) / RAND_MAX;
    if (r > prob) {
      up.push_back(trk);
    }
  }

  return std::make_tuple(central, up);
}

std::pair<std::vector<fastjet::PseudoJet>, std::vector<std::vector<fastjet::PseudoJet>>> getAk15Jets(
    std::vector<RecTrackFormat> tracks,
    std::vector<RecLeptonFormat> leptons,
    double pt_cut,
    double eta_cut,
    double d0_cut,
    double dz_cut,
    double dr_cut)
{
  std::vector<fastjet::PseudoJet> input_particles;

  for (const auto &track : tracks)
  {

    if (track.pt() >= pt_cut &&
        std::abs(track.eta()) <= eta_cut &&
        std::abs(track.d0()) < d0_cut &&
        std::abs(track.dz()) < dz_cut &&
        track.dr(leptons.at(0)) >= dr_cut &&
        track.dr(leptons.at(1)) >= dr_cut)
    {

      double px = track.px();
      double py = track.py();
      double pz = track.pz();
      double p2 = px * px + py * py + pz * pz;
      double mass = getMassFromPDG(track.pdgid());
      double E_corrected = std::sqrt(p2 + mass * mass);
      fastjet::PseudoJet particle(px, py, pz, E_corrected); // Correcting energy of tracks to account for momentum smearing
      input_particles.emplace_back(particle);
    }
  }

  fastjet::JetDefinition jet_def(fastjet::antikt_algorithm, 1.5);
  fastjet::ClusterSequence cluster_seq(input_particles, jet_def);

  std::vector<fastjet::PseudoJet> inclusive_jets = sorted_by_pt(cluster_seq.inclusive_jets(0.0));

  std::vector<std::vector<fastjet::PseudoJet>> cluster_constituents;

  for (const auto &jet : inclusive_jets)
  {
    cluster_constituents.push_back(jet.constituents());
  }

  return std::make_pair(inclusive_jets, cluster_constituents);
}

///////////////////////////////////////////////////////////////
//                       Initialize                          //
// Function called one time at the beginning of the analysis //
///////////////////////////////////////////////////////////////

MAbool user::Initialize(const MA5::Configuration &cfg,
                        const std::map<std::string, std::string> &parameters)
{
  
  // Initializing PhysicsService for RECO
  PHYSICS->recConfig().Reset();

  // ===== Signal region ===== //
  Manager()->AddRegionSelection("SR");
  Manager()->AddRegionSelection("CRDY");
  Manager()->AddRegionSelection("CRTT");

  Manager()->AddRegionSelection("SR_TRACKUP");
  Manager()->AddRegionSelection("CRDY_TRACKUP");
  Manager()->AddRegionSelection("CRTT_TRACKUP");

  // ===== Selections ===== //
  Manager()->AddCut("twoOSSFLeptons");

  Manager()->AddCut("oneAk15Cluster_SR", "SR");
  Manager()->AddCut("oneAk15Cluster_SR_TRACKUP", "SR_TRACKUP");
  Manager()->AddCut("oneAk15Cluster_DY", "CRDY");
  Manager()->AddCut("oneAk15Cluster_DY_TRACKUP", "CRDY_TRACKUP");
  Manager()->AddCut("oneAk15Cluster_TT", "CRTT");
  Manager()->AddCut("oneAk15Cluster_TT_TRACKUP", "CRTT_TRACKUP");

  Manager()->AddCut("onShellZMass_SR", "SR");
  Manager()->AddCut("onShellZMass_SR_TRACKUP", "SR_TRACKUP");
  Manager()->AddCut("offShellZMass_DY", "CRDY");
  Manager()->AddCut("offShellZMass_DY_TRACKUP", "CRDY_TRACKUP");
  Manager()->AddCut("offShellZMass_TT", "CRTT");
  Manager()->AddCut("offShellZMass_TT_TRACKUP", "CRTT_TRACKUP");

  Manager()->AddCut("minZpT");

  Manager()->AddCut("noBTag_SR", "SR");
  Manager()->AddCut("noBTag_SR_TRACKUP", "SR_TRACKUP");
  Manager()->AddCut("noBTag_DY", "CRDY");
  Manager()->AddCut("noBTag_DY_TRACKUP", "CRDY_TRACKUP");
  Manager()->AddCut("yesBTag_TT", "CRTT");
  Manager()->AddCut("yesBTag_TT_TRACKUP", "CRTT_TRACKUP");

  Manager()->AddCut("minAk15pT_SR", "SR");
  Manager()->AddCut("minAk15pT_SR_TRACKUP", "SR_TRACKUP");
  Manager()->AddCut("minAk15pT_DY", "CRDY");
  Manager()->AddCut("minAk15pT_DY_TRACKUP", "CRDY_TRACKUP");
  Manager()->AddCut("minAk15pT_TT", "CRTT");
  Manager()->AddCut("minAk15pT_TT_TRACKUP", "CRTT_TRACKUP");

  Manager()->AddCut("dROverlap_SR", "SR");
  Manager()->AddCut("dROverlap_SR_TRACKUP", "SR_TRACKUP");
  Manager()->AddCut("dROverlap_DY", "CRDY");
  Manager()->AddCut("dROverlap_DY_TRACKUP", "CRDY_TRACKUP");
  Manager()->AddCut("dROverlap_TT", "CRTT");
  Manager()->AddCut("dROverlap_TT_TRACKUP", "CRTT_TRACKUP");

  // ===== Histograms ===== //

  for (const auto& r : regionsForHistos)
  {
    // Z Histograms
    Manager()->AddHisto("ZMass_"+r, 60, 0.0, 300.0, r);
    Manager()->AddHisto("ZpT_"+r, 500, 0.0, 500.0, r);
    Manager()->AddHisto("ZETA_"+r, 40, -5.0, 5.0, r);
    Manager()->AddHisto("ZPHI_"+r, 40, -3.14, 3.14, r);

    // Lepton Histograms
    Manager()->AddHisto("LeadLepPT_"+r, 300, 0.0, 300.0, r);
    Manager()->AddHisto("SubleadLepPT_"+r, 300, 0.0, 300.0, r);
    Manager()->AddHisto("LeadLepETA_"+r, 20, -2.5, 2.5, r);
    Manager()->AddHisto("SubleadLepETA_"+r, 20, -2.5, 2.5, r);
    Manager()->AddHisto("LeadLepPHI_"+r, 20, -3.14, 3.14, r);
    Manager()->AddHisto("SubleadLepPHI_"+r, 20, -3.14, 3.14, r);

    // Muon Histograms
    Manager()->AddHisto("LeadMuonPT_"+r, 300, 0.0, 300.0, r);
    Manager()->AddHisto("SubleadMuonPT_"+r, 300, 0.0, 300.0, r);
    Manager()->AddHisto("LeadMuonETA_"+r, 20, -2.5, 2.5, r);
    Manager()->AddHisto("SubleadMuonETA_"+r, 20, -2.5, 2.5, r);
    Manager()->AddHisto("LeadMuonPHI_"+r, 20, -3.14, 3.14, r);
    Manager()->AddHisto("SubleadMuonPHI_"+r, 20, -3.14, 3.14, r);

    // Electron Histograms
    Manager()->AddHisto("LeadElecPT_"+r, 300, 0.0, 300.0, r);
    Manager()->AddHisto("SubleadElecPT_"+r, 300, 0.0, 300.0, r);
    Manager()->AddHisto("LeadElecETA_"+r, 20, -2.5, 2.5, r);
    Manager()->AddHisto("SubleadElecETA_"+r, 20, -2.5, 2.5, r);
    Manager()->AddHisto("LeadElecPHI_"+r, 20, -3.14, 3.14, r);
    Manager()->AddHisto("SubleadElecPHI_"+r, 20, -3.14, 3.14, r);

    // Ak4 Histograms
    Manager()->AddHisto("NJets_"+r, 10, 0.0, 10.0, r);
    Manager()->AddHisto("LeadAk4PT_"+r, 500, 0.0, 500.0, r);
    Manager()->AddHisto("LeadAk4ETA_"+r, 100, -5.0, 5.0, r);
    Manager()->AddHisto("LeadAk4PHI_"+r, 40, -3.14, 3.14, r);
    Manager()->AddHisto("LeadAk4NTRACKS_"+r, 100, 0.0, 100.0, r);

    // Ak15 Histograms
    Manager()->AddHisto("LeadAk15PT_"+r, 500, 0.0, 500.0, r);
    Manager()->AddHisto("LeadAk15NTRACKS_"+r, 200, 0.0, 200.0, r);
    Manager()->AddHisto("LeadAk15ETA_"+r, 40, -5.0, 5.0, r);
    Manager()->AddHisto("LeadAk15PHI_"+r, 40, -3.14, 3.14, r);
    Manager()->AddHisto("LeadAk15Mass_"+r, 400, 0, 400, r);

    // Sphericity Histograms
    Manager()->AddHisto("labSphericity_"+r, 50, 0.0, 1.0, r);
    Manager()->AddHisto("boostedSphericity_"+r, 50, 0.0, 1.0, r);

    // MET Histograms
    Manager()->AddHisto("MET_PT_"+r, 500, 0.0, 500.0, r);
    Manager()->AddHisto("MET_PHI_"+r, 40, -3.14, 3.14, r);

    // ABCD Histograms
    Manager()->AddHisto("ABCD_A_"+r, 100, 21.0, 121.0, r);
    Manager()->AddHisto("ABCD_B1_"+r, 100, 21.0, 121.0, r);
    Manager()->AddHisto("ABCD_B2_"+r, 100, 21.0, 121.0, r);
    Manager()->AddHisto("ABCD_C1_"+r, 1, 14.0, 21.0, r);
    Manager()->AddHisto("ABCD_C2_"+r, 1, 0.0, 14.0, r);
    Manager()->AddHisto("ABCD_D1_"+r, 1, 14.0, 21.0, r);
    Manager()->AddHisto("ABCD_D2_"+r, 1, 0.0, 14.0, r);
    Manager()->AddHisto("ABCD_E1_"+r, 1, 14.0, 21.0, r);
    Manager()->AddHisto("ABCD_E2_"+r, 1, 0.0, 14.0, r);
  }

  return true;
}

///////////////////////////////////////////////////////////////
//                          Execute                          //
//        Function called each time an event is read         //
///////////////////////////////////////////////////////////////

bool user::Execute(SampleFormat &sample, const EventFormat &event)
{
  
  // Event weight
  double weight = 1.;
  if (!Configuration().IsNoEventWeight() && event.mc() != 0)
  {
    weight = event.mc()->weight();
  }
  Manager()->InitializeForNewEvent(weight);
  if (event.rec() == 0)
  {
    return true;
  }

  // Eflow Collections
  std::vector<RecTrackFormat> eflowTracks = event.rec()->EFlowTracks();

  // Track dropping systematic: get central and variation
  std::vector<RecTrackFormat> centralTracks, TRACKUP;
  std::tie(centralTracks, TRACKUP) = doTracksDropping(event.rec()->tracks());
  
  std::vector<RecParticleFormat> eflowPhotons = event.rec()->EFlowPhotons();
  std::vector<RecParticleFormat> eflowNeutralHadrons = event.rec()->EFlowNeutralHadrons();

  // DEFINING OBJECT CUTS

  // ELECTRONS
  float const ELECTRON_PT_MIN = 15;
  float const ELECTRON_ETA_MAX = 2.5;
  float const ELECTRON_ISO_PT_MIN = 0.5;
  float const ELECTRON_ISO_DR_MAX = 0.3;
  float const ELECTRON_ISO_DR_MIN = 0.01;

  // MUONS
  float const MUON_PT_MIN = 10;
  float const MUON_ETA_MAX = 2.4;
  float const MUON_DZ = 0.1;
  float const MUON_D0 = 0.02;
  float const MUON_ISO_PT_MIN = 0.1;
  float const MUON_ISO_DR_MAX = 0.4;
  float const MUON_ISO_DR_MIN = 0.01;

  // Ak4jets
  float const AK4JET_PT_MIN = 30.0;
  float const AK4JET_ETA_MAX = 2.5;
  float const AK4LEP_DR = 0.4;

  //////////////////////////////////////////////
  //  Applying Base Lepton Object Selections  //
  //////////////////////////////////////////////

  // Electron Collections
  vector<RecLeptonFormat> electrons = filter_electrons(event.rec()->electrons(), ELECTRON_PT_MIN, ELECTRON_ETA_MAX, ELECTRON_ISO_PT_MIN, ELECTRON_ISO_DR_MAX, ELECTRON_ISO_DR_MIN, centralTracks, eflowPhotons, eflowNeutralHadrons);
  vector<RecLeptonFormat> posElectrons = filter_electrons(electrons, ELECTRON_PT_MIN, ELECTRON_ETA_MAX, ELECTRON_ISO_PT_MIN, ELECTRON_ISO_DR_MAX, ELECTRON_ISO_DR_MIN, centralTracks, eflowPhotons, eflowNeutralHadrons, "+");
  vector<RecLeptonFormat> negElectrons = filter_electrons(electrons, ELECTRON_PT_MIN, ELECTRON_ETA_MAX, ELECTRON_ISO_PT_MIN, ELECTRON_ISO_DR_MAX, ELECTRON_ISO_DR_MIN, centralTracks, eflowPhotons, eflowNeutralHadrons, "-");

  // Muon Collections
  vector<RecLeptonFormat> muons = filter_muons(event.rec()->muons(), MUON_PT_MIN, MUON_ETA_MAX, MUON_D0, MUON_DZ, MUON_ISO_PT_MIN, MUON_ISO_DR_MAX, MUON_ISO_DR_MIN, centralTracks, eflowPhotons, eflowNeutralHadrons);
  vector<RecLeptonFormat> posMuons = filter_muons(muons, MUON_PT_MIN, MUON_ETA_MAX, MUON_D0, MUON_DZ, MUON_ISO_PT_MIN, MUON_ISO_DR_MAX, MUON_ISO_DR_MIN, centralTracks, eflowPhotons, eflowNeutralHadrons, "+");
  vector<RecLeptonFormat> negMuons = filter_muons(muons, MUON_PT_MIN, MUON_ETA_MAX, MUON_D0, MUON_DZ, MUON_ISO_PT_MIN, MUON_ISO_DR_MAX, MUON_ISO_DR_MIN, centralTracks, eflowPhotons, eflowNeutralHadrons, "-");

  vector<RecLeptonFormat> leptons = electrons;               // Start with electrons
  leptons.insert(leptons.end(), muons.begin(), muons.end()); // Add muons

  vector<RecLeptonFormat> posLeptons = posElectrons;                     // Start with positive electrons
  posLeptons.insert(posLeptons.end(), posMuons.begin(), posMuons.end()); // Add positive muons

  vector<RecLeptonFormat> negLeptons = negElectrons;                     // Start with negative electrons
  negLeptons.insert(negLeptons.end(), negMuons.begin(), negMuons.end()); // Add negative muons

  // Sort all collections before event selection
  sortLeptonCollections(electrons, posElectrons, negElectrons, muons, posMuons, negMuons, leptons, posLeptons, negLeptons);

  //////////////////////////////////////////////
  //           Event level selections         //
  //////////////////////////////////////////////

  // Event Cut definitions
  float const LEAD_LEPTON_PT = 25;

  float const ZMASS_LOW = 60.0;
  float const ZMASS_HIGH = 120.0;
  float const ZPT_MIN = 25.0;

  float const TRACK_PT_MIN = 1.0;
  float const TRACK_ETA_MAX = 2.5;
  float const TRACK_D0_MAX = 0.05;
  float const TRACK_DZ_MAX = 0.05;
  float const TRACK_DR_MAX = 0.4;

  float const AK15JET_PT_MIN = 60.0;
  float const DR_MAX = 1.5;

  const double DR_MATCH_SOFTLEP = 0.4;

  // Two OSSF Lepton + Lead Lepton pT Selection
  bool twoOSleptons = (posLeptons.size() == 1 && negLeptons.size() == 1);              // Require exactly two opposite-sign leptons (either muons or electrons)
  bool twoSFleptons = (muons.size() == 2 || electrons.size() == 2);                    // Require exactly two muons or exactly two electrons
  bool LeadpTleptons = (leptons.size() > 0) && (leptons.at(0).pt() >= LEAD_LEPTON_PT); // Require the leading lepton pT to be >= 25 GeV
  bool twoOSSFLeptons = twoOSleptons && twoSFleptons && LeadpTleptons;                 // Concatenating cuts
  if (not Manager()->ApplyCut(twoOSSFLeptons, "twoOSSFLeptons")) {
    return true;
  }

  // Do Ak4 clustering
  auto Ak4result = getAk4jets(eflowTracks, eflowPhotons, eflowNeutralHadrons, leptons, AK4LEP_DR);
  std::vector<fastjet::PseudoJet> Ak4jets = applyEnergyScale(Ak4result.first);
  auto ak4Pair = filter_Ak4jetsAndConstituents(Ak4jets, Ak4result.second, AK4JET_PT_MIN, AK4JET_ETA_MAX, AK4LEP_DR, leptons);
  Ak4jets = ak4Pair.first;
  std::vector<std::vector<fastjet::PseudoJet>> Ak4jetConstituents = ak4Pair.second;

  // Do Ak15 clustering
  auto Ak15result = getAk15Jets(centralTracks, leptons, TRACK_PT_MIN, TRACK_ETA_MAX, TRACK_D0_MAX, TRACK_DZ_MAX, TRACK_DR_MAX);
  std::vector<fastjet::PseudoJet> Ak15jets = Ak15result.first;
  std::vector<std::vector<fastjet::PseudoJet>> Ak15jetConstituents = Ak15result.second;

  // Do Ak15 clustering (TRACKUP)
  auto Ak15result_TRACKUP = getAk15Jets(TRACKUP, leptons, TRACK_PT_MIN, TRACK_ETA_MAX, TRACK_D0_MAX, TRACK_DZ_MAX, TRACK_DR_MAX);
  std::vector<fastjet::PseudoJet> Ak15jets_TRACKUP = Ak15result_TRACKUP.first;
  std::vector<std::vector<fastjet::PseudoJet>> Ak15jetConstituents_TRACKUP = Ak15result_TRACKUP.second;

  // One Ak15 cluster
  bool oneAk15Cluster = (Ak15jets.size() > 0);
  bool oneAk15Cluster_TRACKUP = (Ak15jets_TRACKUP.size() > 0);
  Manager()->ApplyCut(oneAk15Cluster, "oneAk15Cluster_SR");
  Manager()->ApplyCut(oneAk15Cluster_TRACKUP, "oneAk15Cluster_SR_TRACKUP");
  Manager()->ApplyCut(oneAk15Cluster, "oneAk15Cluster_DY");
  Manager()->ApplyCut(oneAk15Cluster_TRACKUP, "oneAk15Cluster_DY_TRACKUP");
  Manager()->ApplyCut(oneAk15Cluster, "oneAk15Cluster_TT");
  Manager()->ApplyCut(oneAk15Cluster_TRACKUP, "oneAk15Cluster_TT_TRACKUP");

  // Reconstruct the Z

  // At this point there are exactly one positive and one negative lepton
  RecLeptonFormat posLepton = posLeptons.at(0);
  RecLeptonFormat negLepton = negLeptons.at(0);
  ParticleBaseFormat recoZ;
  recoZ += posLepton.momentum();
  recoZ += negLepton.momentum();

  // ZMass Selection
  bool onShellZMass = (recoZ.m() >= ZMASS_LOW && recoZ.m() <= ZMASS_HIGH);
  Manager()->ApplyCut(onShellZMass, "onShellZMass_SR");
  Manager()->ApplyCut(onShellZMass, "onShellZMass_SR_TRACKUP");

  bool offShellZMass = (recoZ.m() >= ZMASS_HIGH);
  Manager()->ApplyCut(offShellZMass, "offShellZMass_DY");
  Manager()->ApplyCut(offShellZMass, "offShellZMass_DY_TRACKUP");
  Manager()->ApplyCut(offShellZMass, "offShellZMass_TT");
  Manager()->ApplyCut(offShellZMass, "offShellZMass_TT_TRACKUP");

  // ZpT Selection
  bool minZpT = (recoZ.pt() >= ZPT_MIN);
  if (not Manager()->ApplyCut(minZpT, "minZpT"))
    return true;

  // Compute soft‑lepton mistag SF for central and up tracks
  auto Ak4jetMistagSF         = computeSoftLeptonMistagSF(Ak4jets, centralTracks, DR_MATCH_SOFTLEP);
  bool noBTag = BTagVeto(Ak4jets, Ak4jetMistagSF);
  Manager()->ApplyCut(noBTag, "noBTag_SR");
  Manager()->ApplyCut(noBTag, "noBTag_DY");
  Manager()->ApplyCut(!noBTag, "yesBTag_TT");
  Manager()->ApplyCut(noBTag, "noBTag_SR_TRACKUP");
  Manager()->ApplyCut(noBTag, "noBTag_DY_TRACKUP");
  Manager()->ApplyCut(!noBTag, "yesBTag_TT_TRACKUP");

  // Ak15 pT
  bool minAk15pT = (!Ak15jets.empty() && Ak15jets.at(0).pt() > AK15JET_PT_MIN);
  bool minAk15pT_TRACKUP = (!Ak15jets_TRACKUP.empty() && Ak15jets_TRACKUP.at(0).pt() > AK15JET_PT_MIN);
  Manager()->ApplyCut(minAk15pT, "minAk15pT_SR");
  Manager()->ApplyCut(minAk15pT, "minAk15pT_DY");
  Manager()->ApplyCut(minAk15pT, "minAk15pT_TT");
  Manager()->ApplyCut(minAk15pT_TRACKUP, "minAk15pT_SR_TRACKUP");
  Manager()->ApplyCut(minAk15pT_TRACKUP, "minAk15pT_DY_TRACKUP");
  Manager()->ApplyCut(minAk15pT_TRACKUP, "minAk15pT_TT_TRACKUP");

  bool dROverlap = false;
  bool dROverlap_TRACKUP = false;
  if (Ak4jets.empty())
  {
    Manager()->ApplyCut(dROverlap, "dROverlap_SR");
    Manager()->ApplyCut(dROverlap, "dROverlap_DY");
    Manager()->ApplyCut(dROverlap, "dROverlap_TT");
    Manager()->ApplyCut(dROverlap_TRACKUP, "dROverlap_SR_TRACKUP");
    Manager()->ApplyCut(dROverlap_TRACKUP, "dROverlap_DY_TRACKUP");
    Manager()->ApplyCut(dROverlap_TRACKUP, "dROverlap_TT_TRACKUP");
    return true;
  }
  else
  {
    dROverlap = (!Ak15jets.empty() && calculateDeltaR(Ak4jets.at(0), Ak15jets.at(0)) < DR_MAX);
    dROverlap_TRACKUP = (!Ak15jets_TRACKUP.empty() && calculateDeltaR(Ak4jets.at(0), Ak15jets_TRACKUP.at(0)) < DR_MAX);
  }
  
  Manager()->ApplyCut(dROverlap, "dROverlap_SR");
  Manager()->ApplyCut(dROverlap, "dROverlap_DY");
  Manager()->ApplyCut(dROverlap, "dROverlap_TT");
  Manager()->ApplyCut(dROverlap_TRACKUP, "dROverlap_SR_TRACKUP");
  Manager()->ApplyCut(dROverlap_TRACKUP, "dROverlap_DY_TRACKUP");
  Manager()->ApplyCut(dROverlap_TRACKUP, "dROverlap_TT_TRACKUP");

  //////////////////////////////////////////////
  //           Filling Histograms             //
  //////////////////////////////////////////////
  
  // Calculating Lab and Boosted Sphericity
  double labSphericity_CENTRAL = 0.0;
  double boostedSphericity_CENTRAL = 0.0;
  double labSphericity_TRACKUP = 0.0;
  double boostedSphericity_TRACKUP = 0.0;

  if (!Ak15jets.empty()) {
    labSphericity_CENTRAL = sphericity(Ak15jetConstituents.at(0), 1.0);
    if (Ak15jetConstituents.at(0).size() > 1) {
      auto boosted_CENTRAL = boostToSUEP(Ak15jetConstituents.at(0), Ak15jets.at(0));
      boostedSphericity_CENTRAL = sphericity(boosted_CENTRAL, 1.0);
    } 
    else {
      boostedSphericity_CENTRAL = 0.0;
    }
  }

  if (!Ak15jets_TRACKUP.empty()) {
    labSphericity_TRACKUP = sphericity(Ak15jetConstituents_TRACKUP.at(0), 1.0);
    if (Ak15jetConstituents_TRACKUP.at(0).size() > 1) {
      auto boosted_TRACKUP = boostToSUEP(Ak15jetConstituents_TRACKUP.at(0), Ak15jets_TRACKUP.at(0));
      boostedSphericity_TRACKUP = sphericity(boosted_TRACKUP, 1.0);
    } 
    else {
      boostedSphericity_TRACKUP = 0.0;
    }
  }

  // Z Histograms
  FillHistoAllRegions("ZMass",   recoZ.m(), regionsForHistos);
  FillHistoAllRegions("ZpT",     recoZ.pt(), regionsForHistos);
  FillHistoAllRegions("ZETA",    recoZ.eta(), regionsForHistos);
  FillHistoAllRegions("ZPHI",    recoZ.phi(), regionsForHistos);

  // Lepton Histograms
  FillHistoAllRegions("LeadLepPT", leptons.at(0).pt(), regionsForHistos);
  FillHistoAllRegions("SubleadLepPT", leptons.at(1).pt(), regionsForHistos);
  FillHistoAllRegions("LeadLepETA", leptons.at(0).eta(), regionsForHistos);
  FillHistoAllRegions("SubleadLepETA", leptons.at(1).eta(), regionsForHistos);
  FillHistoAllRegions("LeadLepPHI", leptons.at(0).phi(), regionsForHistos);
  FillHistoAllRegions("SubleadLepPHI", leptons.at(1).phi(), regionsForHistos);

  // Muon or Electron Histograms
  if (muons.size() > 0)
  {
    FillHistoAllRegions("LeadMuonPT", muons.at(0).pt(), regionsForHistos);
    FillHistoAllRegions("SubleadMuonPT", muons.at(1).pt(), regionsForHistos);
    FillHistoAllRegions("LeadMuonETA", muons.at(0).eta(), regionsForHistos);
    FillHistoAllRegions("SubleadMuonETA", muons.at(1).eta(), regionsForHistos);
    FillHistoAllRegions("LeadMuonPHI", muons.at(0).phi(), regionsForHistos);
    FillHistoAllRegions("SubleadMuonPHI", muons.at(1).phi(), regionsForHistos);
  }
  else
  {
    FillHistoAllRegions("LeadElecPT", electrons.at(0).pt(), regionsForHistos);
    FillHistoAllRegions("SubleadElecPT", electrons.at(1).pt(), regionsForHistos);
    FillHistoAllRegions("LeadElecETA", electrons.at(0).eta(), regionsForHistos);
    FillHistoAllRegions("SubleadElecETA", electrons.at(1).eta(), regionsForHistos);
    FillHistoAllRegions("LeadElecPHI", electrons.at(0).phi(), regionsForHistos);
    FillHistoAllRegions("SubleadElecPHI", electrons.at(1).phi(), regionsForHistos);
  }

  // Ak4 Histograms
  FillHistoAllRegions("NJets", Ak4jets.size(), regionsForHistos);
  FillHistoAllRegions("LeadAk4PT", Ak4jets.at(0).pt(), regionsForHistos);
  FillHistoAllRegions("LeadAk4ETA", Ak4jets.at(0).eta(), regionsForHistos);
  FillHistoAllRegions("LeadAk4PHI", normalizePhi(Ak4jets.at(0).phi()), regionsForHistos);
  FillHistoAllRegions("LeadAk4NTRACKS", Ak4jetConstituents.at(0).size(), regionsForHistos);

  // MET Histograms
  FillHistoAllRegions("MET_PT", event.rec()->MET().pt(), regionsForHistos);
  FillHistoAllRegions("MET_PHI", event.rec()->MET().phi(), regionsForHistos);

  // Ak15 Histograms
  if(!Ak15jets.empty()){
    FillHistoAllRegions("LeadAk15PT", Ak15jets.at(0).pt(), regionsForHistos_CENTRAL);
    FillHistoAllRegions("LeadAk15NTRACKS", Ak15jetConstituents.at(0).size(), regionsForHistos_CENTRAL);
    FillHistoAllRegions("LeadAk15ETA", Ak15jets.at(0).eta(), regionsForHistos_CENTRAL);
    FillHistoAllRegions("LeadAk15PHI", normalizePhi(Ak15jets.at(0).phi()), regionsForHistos_CENTRAL);
    FillHistoAllRegions("LeadAk15Mass", Ak15jets.at(0).m(), regionsForHistos_CENTRAL);
    FillHistoAllRegions("labSphericity", labSphericity_CENTRAL, regionsForHistos_CENTRAL);
    FillHistoAllRegions("boostedSphericity", boostedSphericity_CENTRAL, regionsForHistos_CENTRAL);

    // Extended ABCD regions (CENTRAL)
    if ((Ak15jetConstituents.at(0).size() >= 21.0) && (Ak4jets.at(0).pt() <= 135.0))
    {
      FillHistoAllRegions("ABCD_A", Ak15jetConstituents.at(0).size(), regionsForHistos_CENTRAL);
    }
    if ((Ak15jetConstituents.at(0).size() >= 14.0) && (Ak15jetConstituents.at(0).size() < 21.0) && (Ak4jets.at(0).pt() <= 135.0))
    {
      FillHistoAllRegions("ABCD_C1", Ak15jetConstituents.at(0).size(), regionsForHistos_CENTRAL);
    }
    if ((Ak15jetConstituents.at(0).size() < 14.0) && (Ak4jets.at(0).pt() <= 135.0))
    {
      FillHistoAllRegions("ABCD_C2", Ak15jetConstituents.at(0).size(), regionsForHistos_CENTRAL);
    }

    if ((Ak15jetConstituents.at(0).size() >= 21.0) && (Ak4jets.at(0).pt() <= 220.0) && (Ak4jets.at(0).pt() > 135.0))
    {
      FillHistoAllRegions("ABCD_B1", Ak15jetConstituents.at(0).size(), regionsForHistos_CENTRAL);
    }
    if ((Ak15jetConstituents.at(0).size() >= 14.0) && (Ak15jetConstituents.at(0).size() < 21.0) && (Ak4jets.at(0).pt() <= 220.0) && (Ak4jets.at(0).pt() > 135.0))
    {
      FillHistoAllRegions("ABCD_D1", Ak15jetConstituents.at(0).size(), regionsForHistos_CENTRAL);
    }
    if ((Ak15jetConstituents.at(0).size() < 14.0) && (Ak4jets.at(0).pt() <= 220.0) && (Ak4jets.at(0).pt() > 135.0))
    {
      FillHistoAllRegions("ABCD_D2", Ak15jetConstituents.at(0).size(), regionsForHistos_CENTRAL);
    }

    if ((Ak15jetConstituents.at(0).size() >= 21.0) && (Ak4jets.at(0).pt() > 220.0))
    {
      FillHistoAllRegions("ABCD_B2", Ak15jetConstituents.at(0).size(), regionsForHistos_CENTRAL);
    }
    if ((Ak15jetConstituents.at(0).size() >= 14.0) && (Ak15jetConstituents.at(0).size() < 21.0) && (Ak4jets.at(0).pt() > 220.0))
    {
      FillHistoAllRegions("ABCD_E1", Ak15jetConstituents.at(0).size(), regionsForHistos_CENTRAL);
    }
    if ((Ak15jetConstituents.at(0).size() < 14.0) && (Ak4jets.at(0).pt() > 220.0))
    {
      FillHistoAllRegions("ABCD_E2", Ak15jetConstituents.at(0).size(), regionsForHistos_CENTRAL);
    }
  }

  // Ak15 Histograms (TRACKUP)
  if(!Ak15jets_TRACKUP.empty()){
    FillHistoAllRegions("LeadAk15PT", Ak15jets_TRACKUP.at(0).pt(), regionsForHistos_TRACKUP);
    FillHistoAllRegions("LeadAk15NTRACKS", Ak15jetConstituents_TRACKUP.at(0).size(), regionsForHistos_TRACKUP);
    FillHistoAllRegions("LeadAk15ETA", Ak15jets_TRACKUP.at(0).eta(), regionsForHistos_TRACKUP);
    FillHistoAllRegions("LeadAk15PHI", normalizePhi(Ak15jets_TRACKUP.at(0).phi()), regionsForHistos_TRACKUP);
    FillHistoAllRegions("LeadAk15Mass", Ak15jets_TRACKUP.at(0).m(), regionsForHistos_TRACKUP);
    FillHistoAllRegions("labSphericity", labSphericity_TRACKUP, regionsForHistos_TRACKUP);
    FillHistoAllRegions("boostedSphericity", boostedSphericity_TRACKUP, regionsForHistos_TRACKUP);

    // Extended ABCD regions (TRACKUP)
    if ((Ak15jetConstituents_TRACKUP.at(0).size() >= 21.0) && (Ak4jets.at(0).pt() <= 135.0))
    {
      FillHistoAllRegions("ABCD_A", Ak15jetConstituents_TRACKUP.at(0).size(), regionsForHistos_TRACKUP);
    }
    if ((Ak15jetConstituents_TRACKUP.at(0).size() >= 14.0) && (Ak15jetConstituents_TRACKUP.at(0).size() < 21.0) && (Ak4jets.at(0).pt() <= 135.0))
    {
      FillHistoAllRegions("ABCD_C1", Ak15jetConstituents_TRACKUP.at(0).size(), regionsForHistos_TRACKUP);
    }
    if ((Ak15jetConstituents_TRACKUP.at(0).size() < 14.0) && (Ak4jets.at(0).pt() <= 135.0))
    {
      FillHistoAllRegions("ABCD_C2", Ak15jetConstituents_TRACKUP.at(0).size(), regionsForHistos_TRACKUP);
    }

    if ((Ak15jetConstituents_TRACKUP.at(0).size() >= 21.0) && (Ak4jets.at(0).pt() <= 220.0) && (Ak4jets.at(0).pt() > 135.0))
    {
      FillHistoAllRegions("ABCD_B1", Ak15jetConstituents_TRACKUP.at(0).size(), regionsForHistos_TRACKUP);
    }
    if ((Ak15jetConstituents_TRACKUP.at(0).size() >= 14.0) && (Ak15jetConstituents_TRACKUP.at(0).size() < 21.0) && (Ak4jets.at(0).pt() <= 220.0) && (Ak4jets.at(0).pt() > 135.0))
    {
      FillHistoAllRegions("ABCD_D1", Ak15jetConstituents_TRACKUP.at(0).size(), regionsForHistos_TRACKUP);
    }
    if ((Ak15jetConstituents_TRACKUP.at(0).size() < 14.0) && (Ak4jets.at(0).pt() <= 220.0) && (Ak4jets.at(0).pt() > 135.0))
    {
      FillHistoAllRegions("ABCD_D2", Ak15jetConstituents_TRACKUP.at(0).size(), regionsForHistos_TRACKUP);
    }

    if ((Ak15jetConstituents_TRACKUP.at(0).size() >= 21.0) && (Ak4jets.at(0).pt() > 220.0))
    {
      FillHistoAllRegions("ABCD_B2", Ak15jetConstituents_TRACKUP.at(0).size(), regionsForHistos_TRACKUP);
    }
    if ((Ak15jetConstituents_TRACKUP.at(0).size() >= 14.0) && (Ak15jetConstituents_TRACKUP.at(0).size() < 21.0) && (Ak4jets.at(0).pt() > 220.0))
    {
      FillHistoAllRegions("ABCD_E1", Ak15jetConstituents_TRACKUP.at(0).size(), regionsForHistos_TRACKUP);
    }
    if ((Ak15jetConstituents_TRACKUP.at(0).size() < 14.0) && (Ak4jets.at(0).pt() > 220.0))
    {
      FillHistoAllRegions("ABCD_E2", Ak15jetConstituents_TRACKUP.at(0).size(), regionsForHistos_TRACKUP);
    }

  }

  return true;
}

///////////////////////////////////////////////////////////////
//                        Finalize                           //
//    Function called one time at the end of the analysis    //
///////////////////////////////////////////////////////////////
void user::Finalize(const SampleFormat &summary, const std::vector<SampleFormat> &files)
{
}

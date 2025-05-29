#include "SampleAnalyzer/User/Analyzer/user.h"
using namespace MA5;
using namespace std;
#include <Eigen/Dense>
#include <random>
#include <cmath>

const std::vector<std::string> regionsForHistos = {"SR", "SR_TRACKUP"};
const std::vector<std::string> regionsForHistos_CENTRAL = {"SR"};
const std::vector<std::string> regionsForHistos_TRACKUP = {"SR_TRACKUP"};

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

template <typename T, typename U>
bool dPhiCut(const T &obj1, const U &obj2, float minDphi)
{
  float dphi = computeDeltaPhi(obj1.phi(), obj2.phi());
  return std::fabs(dphi) > minDphi;
}

template <typename T, typename U>
double dPhi_Ak4_MET(const T &met, const std::vector<U> &ak4jets)
{
  std::vector<double> dPhiValues;
  for (const auto &jet : ak4jets)
  {
    float dphi = computeDeltaPhi(met.phi(), jet.phi());
    dPhiValues.push_back(std::fabs(dphi));
  }
  std::sort(dPhiValues.begin(), dPhiValues.end());
  double minDphiValue = dPhiValues.front();

  return minDphiValue;
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
    throw std::runtime_error("[FATAL] Unrecognized PDG id in getMassFromPDG: " + std::to_string(pdgid) + " - please update the function with this particle's mass.");
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
    {
      totalpt += pt;
    }
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

vector<RecLeptonFormat> filter_muons(const vector<RecLeptonFormat> &objects,
                                     float ptmin,
                                     float etamax,
                                     float d0,
                                     float dz,
                                     float iso_pTMin,
                                     float iso_dRMax,
                                     float iso_dRMin,
                                     const vector<RecTrackFormat> &tracks,
                                     const vector<RecParticleFormat> &eflowPhotons,
                                     const vector<RecParticleFormat> &eflowNeutralHadrons,
                                     const string &selection,
                                     string charge = "")
{

  // Helper function to select muons
  vector<RecLeptonFormat> filtered;
  for (const auto &obj : objects)
  {
    // Charge selection
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

    // Independent of loose or tight
    if (fabs(obj.eta()) > etamax)
      continue;

    if (fabs(obj.d0()) > d0)
      continue;

    // Loose or tight dependence
    if (obj.pt() < ptmin)
      continue;

    if (fabs(obj.dz()) > dz)
      continue;

    // Iso selection: If loose, use nominal Delphes output. If tight, calculate on the fly
    if (selection == "loose")
    {
      if (!iso(obj, tracks, eflowPhotons, eflowNeutralHadrons, iso_pTMin, iso_dRMax, iso_dRMin, "pfIso2"))
        continue;
      else
      {
        filtered.push_back(obj);
      }
    }
    else if (selection == "tight")
    {
      if (!iso(obj, tracks, eflowPhotons, eflowNeutralHadrons, iso_pTMin, iso_dRMax, iso_dRMin, "pfIso5"))
        continue;
      else
      {
        filtered.push_back(obj);
      }
    }
  }

  return filtered;
}

vector<RecLeptonFormat> filter_electrons(const vector<RecLeptonFormat> &objects,
                                         float ptmin,
                                         float etamax,
                                         float iso_pTMin,
                                         float iso_dRMax,
                                         float iso_dRMin,
                                         const vector<RecTrackFormat> &centralTracks,
                                         const vector<RecParticleFormat> &eflowPhotons,
                                         const vector<RecParticleFormat> &eflowNeutralHadrons,
                                         const string &selection,
                                         string charge = "")
{
  // Helper function to select electrons
  vector<RecLeptonFormat> filtered;

  for (const auto &obj : objects)
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

    // Loose or tight dependence
    if (obj.pt() < ptmin)
      continue;

    // Independent of loose or right
    if (fabs(obj.eta()) > etamax)
      continue;

    // Exclude crack region: electrons with 1.444 < |eta| < 1.566
    if ((fabs(obj.eta()) > 1.444) && (fabs(obj.eta()) < 1.566))
      continue;

    if (fabs(obj.d0()) > (0.05 + 0.05 * (fabs(obj.eta()) > 1.479)))
      continue;

    if (fabs(obj.dz()) > (0.1 + 0.1 * (fabs(obj.eta()) > 1.479)))
      continue;

    // Iso selection: If loose, use nominal Delphes output. If tight, calculate on the fly
    if (selection == "loose")
    {
      if (!iso(obj, centralTracks, eflowPhotons, eflowNeutralHadrons, iso_pTMin, iso_dRMax, iso_dRMin, "WP90"))
        continue;
      else
      {
        filtered.push_back(obj);
      }
    }
    else if (selection == "tight")
    {
      if (!iso(obj, centralTracks, eflowPhotons, eflowNeutralHadrons, iso_pTMin, iso_dRMax, iso_dRMin, "WP80"))
        continue;
      else
      {
        filtered.push_back(obj);
      }
    }
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

  for (std::size_t i = 0; i < jets.size(); ++i)
  {
    const auto &jet = jets[i];
    double pt = jet.pt();
    double eff = 0.01 + 0.000038 * pt;

    if (i < mistagSF.size())
    {
      eff = std::min(1.0, eff * mistagSF[i]);
    }

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
    // ------------------------------------
    // probabilistic track‑inefficiency test
    // ------------------------------------
    //int absId = std::abs(track.pdgid());
    //double dropProb = (absId == 11 || absId == 13) ? 0.06 : 0.02;
  
    //double rnd = static_cast<double>(rand()) / RAND_MAX;
    //if (rnd < dropProb) continue;        // skip this track
  
    // ------------------------------------
    // nominal quality cuts
    // ------------------------------------
    if (track.pt()  <  pt_cut)              continue;
    if (std::abs(track.eta()) > eta_cut)    continue;
    if (std::abs(track.d0())  >= d0_cut)    continue;
    if (std::abs(track.dz())  >= dz_cut)    continue;
    if (track.dr(leptons.at(0)) < dr_cut)   continue;
  
    // ------------------------------------
    // build the PseudoJet and store
    // ------------------------------------
    double px = track.px();
    double py = track.py();
    double pz = track.pz();
    double p2 = px*px + py*py + pz*pz;
    double mass = getMassFromPDG(track.pdgid());
    double E_corrected = std::sqrt(p2 + mass*mass);
  
    fastjet::PseudoJet particle(px, py, pz, E_corrected);
    input_particles.emplace_back(particle);
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
  Manager()->AddRegionSelection("SR_TRACKUP");

  // ===== Selections ===== //
  Manager()->AddCut("orthogonality");
  Manager()->AddCut("oneTightLep");

  Manager()->AddCut("Ak15", "SR");
  Manager()->AddCut("Ak15_TRACKUP", "SR_TRACKUP");

  Manager()->AddCut("oneAk4");

  // ===== Histograms ===== //

  for (const auto& r : regionsForHistos)
  {
    // W Histograms
    Manager()->AddHisto("wTransverseMass_"+r, 100, 30.0, 130.0, r);
    Manager()->AddHisto("wPt_"+r, 94, 60, 1000, r);
    Manager()->AddHisto("wPhi_"+r, 60, -3.2, 3.2, r);

    // Lepton Histograms
    Manager()->AddHisto("lepPt_"+r, 370, 30, 400, r);
    Manager()->AddHisto("lepEta_"+r, 50, -2.5, 2.5, r);
    Manager()->AddHisto("lepPhi_"+r, 60, -3.2, 3.2, r);
    Manager()->AddHisto("looseLep_"+r, 10, 0, 10, r);
    Manager()->AddHisto("tightLep_"+r, 10, 0, 10, r);

    // Muon Histograms
    Manager()->AddHisto("muPt_"+r, 370, 30, 400, r);
    Manager()->AddHisto("muEta_"+r, 50, -2.5, 2.5, r);
    Manager()->AddHisto("muPhi_"+r, 60, -3.2, 3.2, r);
    Manager()->AddHisto("looseMu_"+r, 10, 0, 10, r);
    Manager()->AddHisto("tightMu_"+r, 10, 0, 10, r);

    // Electron Histograms
    Manager()->AddHisto("elePt_"+r, 370, 30, 400, r);
    Manager()->AddHisto("eleEta_"+r, 50, -2.5, 2.5, r);
    Manager()->AddHisto("elePhi_"+r, 60, -3.2, 3.2, r);
    Manager()->AddHisto("looseEle_"+r, 10, 0, 10, r);
    Manager()->AddHisto("tightEle_"+r, 10, 0, 10, r);

    // Ak4 Histograms
    Manager()->AddHisto("NJets_"+r, 20, 0.0, 20.0, r);
    Manager()->AddHisto("ak4Pt_"+r, 300, 0, 300, r);
    Manager()->AddHisto("ak4Eta_"+r, 100, -5.0, 5.0, r);
    Manager()->AddHisto("ak4Phi_"+r, 60, -3.2, 3.2, r);
    Manager()->AddHisto("ak4NTracks_"+r, 100, 0.0, 100.0, r);

    // Ak15 Histograms
    Manager()->AddHisto("ak15Pt_"+r, 500, 0, 500, r);
    Manager()->AddHisto("ak15Eta_"+r, 50, -2.5, 2.5, r);
    Manager()->AddHisto("ak15Phi_"+r, 50, -3.25, 3.25, r);
    Manager()->AddHisto("ak15NTracks_"+r, 200, 0.0, 200.0, r);
    Manager()->AddHisto("ak15Mass_"+r, 30, 0, 400, r);

    // Sphericity Histograms
    Manager()->AddHisto("labSphericity_"+r, 70, 0.3, 1.0, r);
    Manager()->AddHisto("boostedSphericity_"+r, 70, 0.3, 1.0, r);

    // MET Histograms
    Manager()->AddHisto("metPt_"+r, 94, 30, 500, r);
    Manager()->AddHisto("metPhi_"+r, 100, -3, 3, r);

    Manager()->AddHisto("dPhi_MET_Lep_"+r, 100, -3.14, 3.14, r);
    Manager()->AddHisto("dPhi_MET_SUEP_"+r, 100, -3.14, 3.14, r);

    // ABCD Histograms
    Manager()->AddHisto("ABCD_A_"+r, 10, 10.0, 20.0, r);
    Manager()->AddHisto("ABCD_B_"+r, 10, 20.0, 30.0, r);
    Manager()->AddHisto("ABCD_C_"+r, 170, 30.0, 200.0, r);
    Manager()->AddHisto("ABCD_D_"+r, 10, 10.0, 20.0, r);
    Manager()->AddHisto("ABCD_E_"+r, 10, 20.0, 30.0, r);
    Manager()->AddHisto("ABCD_F0_"+r, 10, 30.0, 40.0, r);
    Manager()->AddHisto("ABCD_F1_"+r, 10, 40.0, 50.0, r);
    Manager()->AddHisto("ABCD_F2_"+r, 10, 50.0, 60.0, r);
    Manager()->AddHisto("ABCD_F3_"+r, 20, 60.0, 80.0, r);
    Manager()->AddHisto("ABCD_F4_"+r, 120, 80.0, 200.0, r);
    Manager()->AddHisto("ABCD_G_"+r, 10, 10.0, 20.0, r);
    Manager()->AddHisto("ABCD_H_"+r, 10, 20.0, 30.0, r);
    Manager()->AddHisto("ABCD_SR0_"+r, 10, 30.0, 40.0, r);
    Manager()->AddHisto("ABCD_SR1_"+r, 10, 40.0, 50.0, r);
    Manager()->AddHisto("ABCD_SR2_"+r, 10, 50.0, 60.0, r);
    Manager()->AddHisto("ABCD_SR3_"+r, 20, 60.0, 80.0, r);
    Manager()->AddHisto("ABCD_SR4_"+r, 120, 80.0, 200.0, r);
  
    // ABCD Histograms
    if(r == "SR"){
        Manager()->AddHisto("ABCD_A_ZH_"+r, 100, 21.0, 121.0, r);
        Manager()->AddHisto("ABCD_B1_ZH_"+r, 100, 21.0, 121.0, r);
        Manager()->AddHisto("ABCD_B2_ZH_"+r, 100, 21.0, 121.0, r);
        Manager()->AddHisto("ABCD_C1_ZH_"+r, 1, 14.0, 21.0, r);
        Manager()->AddHisto("ABCD_C2_ZH_"+r, 1, 0.0, 14.0, r);
        Manager()->AddHisto("ABCD_D1_ZH_"+r, 1, 14.0, 21.0, r);
        Manager()->AddHisto("ABCD_D2_ZH_"+r, 1, 0.0, 14.0, r);
        Manager()->AddHisto("ABCD_E1_ZH_"+r, 1, 14.0, 21.0, r);
        Manager()->AddHisto("ABCD_E2_ZH_"+r, 1, 0.0, 14.0, r);
    }

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
  
  // Copy REC MET and apply positive-only fractional smearing
  RecParticleFormat smearedMET = event.rec()->MET();
  double orig_pt = smearedMET.pt();
  double phi    = smearedMET.phi();
  // Random-number generator (thread-local)
  static thread_local std::mt19937 gen(std::random_device{}());
  // Compute dynamic μ and σ with your provided coefficients
  double mu    = 21.334   / (46.505   + 0.53795 * orig_pt);
  double sigma = 3.954    / (5.340    + 0.07439 * orig_pt);
  // Draw fractional smear from N(μ, σ)
  std::normal_distribution<double> dist_rel(mu, sigma);
  double delta_rel = dist_rel(gen);
  // Apply smear only if positive
  double pt_smeared  = (delta_rel > 0.0) ? orig_pt * (1.0 + delta_rel) : orig_pt;
  double px_smeared  = pt_smeared * std::cos(phi);
  double py_smeared  = pt_smeared * std::sin(phi);
  double E_smeared   = pt_smeared;
  smearedMET.momentum().SetPxPyPzE(px_smeared, py_smeared, 0.0, E_smeared);

  // DEFINING OBJECT CUTS

  // Isolation Eflow Collections
  std::vector<RecTrackFormat> eflowTracks = event.rec()->EFlowTracks();

  // Track dropping systematic: get central and variation
  std::vector<RecTrackFormat> centralTracks, TRACKUP;
  std::tie(centralTracks, TRACKUP) = doTracksDropping(event.rec()->tracks());

  std::vector<RecParticleFormat> eflowPhotons = event.rec()->EFlowPhotons();
  std::vector<RecParticleFormat> eflowNeutralHadrons = event.rec()->EFlowNeutralHadrons();

  // ELECTRONS
  float const LOOSE_ELECTRON_PT_MIN = 15;
  float const TIGHT_ELECTRON_PT_MIN = 35;
  float const ELECTRON_ETA_MAX = 2.5;
  float const ELECTRON_ISO_PT_MIN = 0.5;
  float const ELECTRON_ISO_DR_MAX = 0.3;
  float const ELECTRON_ISO_DR_MIN = 0.01;

  // MUONS
  float const LOOSE_MUON_PT_MIN = 10;
  float const TIGHT_MUON_PT_MIN = 30;
  float const MUON_ETA_MAX = 2.4;
  float const LOOSE_MUON_DZ = 0.1;
  float const TIGHT_MUON_DZ = 0.05;
  float const MUON_D0 = 0.02;
  float const MUON_ISO_PT_MIN = 0.1;
  float const MUON_ISO_DR_MAX = 0.4;
  float const MUON_ISO_DR_MIN = 0.01;

  // Ak4jets
  float const AK4_PT_MIN = 30;
  float const AK4_ETA_MAX = 2.4;
  float const AK4LEP_DR = 0.4;
  const double DR_MATCH_SOFTLEP = 0.4;

  //////////////////////////////////////////////
  //  Applying Base Lepton Object Selections  //
  //////////////////////////////////////////////

  // Electron Collections
  vector<RecLeptonFormat> loose_electrons = filter_electrons(event.rec()->electrons(), LOOSE_ELECTRON_PT_MIN, ELECTRON_ETA_MAX, ELECTRON_ISO_PT_MIN, ELECTRON_ISO_DR_MAX, ELECTRON_ISO_DR_MIN, centralTracks, eflowPhotons, eflowNeutralHadrons, "loose");
  vector<RecLeptonFormat> loose_posElectrons = filter_electrons(loose_electrons, LOOSE_ELECTRON_PT_MIN, ELECTRON_ETA_MAX, ELECTRON_ISO_PT_MIN, ELECTRON_ISO_DR_MAX, ELECTRON_ISO_DR_MIN, centralTracks, eflowPhotons, eflowNeutralHadrons, "loose", "+");
  vector<RecLeptonFormat> loose_negElectrons = filter_electrons(loose_electrons, LOOSE_ELECTRON_PT_MIN, ELECTRON_ETA_MAX, ELECTRON_ISO_PT_MIN, ELECTRON_ISO_DR_MAX, ELECTRON_ISO_DR_MIN, centralTracks, eflowPhotons, eflowNeutralHadrons, "loose", "-");

  vector<RecLeptonFormat> tight_electrons = filter_electrons(loose_electrons, TIGHT_ELECTRON_PT_MIN, ELECTRON_ETA_MAX, ELECTRON_ISO_PT_MIN, ELECTRON_ISO_DR_MAX, ELECTRON_ISO_DR_MIN, centralTracks, eflowPhotons, eflowNeutralHadrons, "tight");
  vector<RecLeptonFormat> tight_posElectrons = filter_electrons(tight_electrons, TIGHT_ELECTRON_PT_MIN, ELECTRON_ETA_MAX, ELECTRON_ISO_PT_MIN, ELECTRON_ISO_DR_MAX, ELECTRON_ISO_DR_MIN, centralTracks, eflowPhotons, eflowNeutralHadrons, "tight", "+");
  vector<RecLeptonFormat> tight_negElectrons = filter_electrons(tight_electrons, TIGHT_ELECTRON_PT_MIN, ELECTRON_ETA_MAX, ELECTRON_ISO_PT_MIN, ELECTRON_ISO_DR_MAX, ELECTRON_ISO_DR_MIN, centralTracks, eflowPhotons, eflowNeutralHadrons, "tight", "-");

  // Muon Collections
  vector<RecLeptonFormat> loose_muons = filter_muons(event.rec()->muons(), LOOSE_MUON_PT_MIN, MUON_ETA_MAX, MUON_D0, LOOSE_MUON_DZ, MUON_ISO_PT_MIN, MUON_ISO_DR_MAX, MUON_ISO_DR_MIN, centralTracks, eflowPhotons, eflowNeutralHadrons, "loose");
  vector<RecLeptonFormat> loose_posMuons = filter_muons(loose_muons, LOOSE_MUON_PT_MIN, MUON_ETA_MAX, MUON_D0, LOOSE_MUON_DZ, MUON_ISO_PT_MIN, MUON_ISO_DR_MAX, MUON_ISO_DR_MIN, centralTracks, eflowPhotons, eflowNeutralHadrons, "loose", "+");
  vector<RecLeptonFormat> loose_negMuons = filter_muons(loose_muons, LOOSE_MUON_PT_MIN, MUON_ETA_MAX, MUON_D0, LOOSE_MUON_DZ, MUON_ISO_PT_MIN, MUON_ISO_DR_MAX, MUON_ISO_DR_MIN, centralTracks, eflowPhotons, eflowNeutralHadrons, "loose", "-");

  vector<RecLeptonFormat> tight_muons = filter_muons(loose_muons, TIGHT_MUON_PT_MIN, MUON_ETA_MAX, MUON_D0, TIGHT_MUON_DZ, MUON_ISO_PT_MIN, MUON_ISO_DR_MAX, MUON_ISO_DR_MIN, centralTracks, eflowPhotons, eflowNeutralHadrons, "tight");
  vector<RecLeptonFormat> tight_posMuons = filter_muons(tight_muons, TIGHT_MUON_PT_MIN, MUON_ETA_MAX, MUON_D0, TIGHT_MUON_DZ, MUON_ISO_PT_MIN, MUON_ISO_DR_MAX, MUON_ISO_DR_MIN, centralTracks, eflowPhotons, eflowNeutralHadrons, "tight", "+");
  vector<RecLeptonFormat> tight_negMuons = filter_muons(tight_muons, TIGHT_MUON_PT_MIN, MUON_ETA_MAX, MUON_D0, TIGHT_MUON_DZ, MUON_ISO_PT_MIN, MUON_ISO_DR_MAX, MUON_ISO_DR_MIN, centralTracks, eflowPhotons, eflowNeutralHadrons, "tight", "-");

  // Combining loose Leptons
  vector<RecLeptonFormat> loose_leptons = loose_electrons;
  loose_leptons.insert(loose_leptons.end(), loose_muons.begin(), loose_muons.end());

  vector<RecLeptonFormat> loose_posLeptons = loose_posElectrons;
  loose_posLeptons.insert(loose_posLeptons.end(), loose_posMuons.begin(), loose_posMuons.end());

  vector<RecLeptonFormat> loose_negLeptons = loose_negElectrons;
  loose_negLeptons.insert(loose_negLeptons.end(), loose_negMuons.begin(), loose_negMuons.end());

  // Combining tight Leptons
  vector<RecLeptonFormat> tight_leptons = tight_electrons;
  tight_leptons.insert(tight_leptons.end(), tight_muons.begin(), tight_muons.end());

  vector<RecLeptonFormat> tight_posLeptons = tight_posElectrons;
  tight_posLeptons.insert(tight_posLeptons.end(), tight_posMuons.begin(), tight_posMuons.end());

  vector<RecLeptonFormat> tight_negLeptons = tight_negElectrons;
  tight_negLeptons.insert(tight_negLeptons.end(), tight_negMuons.begin(), tight_negMuons.end());

  // Sorting leptons
  sortLeptonCollections(loose_electrons, loose_posElectrons, loose_negElectrons, loose_muons, loose_posMuons, loose_negMuons, loose_leptons, loose_posLeptons, loose_negLeptons);
  sortLeptonCollections(tight_electrons, tight_posElectrons, tight_negElectrons, tight_muons, tight_posMuons, tight_negMuons, tight_leptons, tight_posLeptons, tight_negLeptons);

  //////////////////////////////////////////////
  //   First Round of Event level selections  //
  //////////////////////////////////////////////

  // Event Cut definitions

  // Orthogonality
  float const ORTHOG_LEAD_LEPTON_PT = 25;

  // MET
  float const FIRST_MET_PT = 20;
  float const SECOND_MET_PT = 30;

  // Orthogonality to GGF offline: Remove events with no leptons
  bool atLeastOneEle = (loose_electrons.size() > 0) && (loose_electrons.at(0).pt() >= ORTHOG_LEAD_LEPTON_PT);
  bool atLeastOneMu = (loose_muons.size() > 0) && (loose_muons.at(0).pt() >= ORTHOG_LEAD_LEPTON_PT);
  bool GGFOrthogonality = (atLeastOneEle || atLeastOneMu);

  // Orthogonality to ZH: Remove events with a pair of loose OSSF leptons
  bool twoOSleptons = (loose_posLeptons.size() == 1 && loose_negLeptons.size() == 1);                     // Require exactly two opposite-sign leptons (either muons or electrons)
  bool twoSFleptons = (loose_muons.size() == 2 || loose_electrons.size() == 2);                           // Require exactly two muons or exactly two electrons
  bool LeadpTleptons = (loose_leptons.size() > 0) && (loose_leptons.at(0).pt() >= ORTHOG_LEAD_LEPTON_PT); // Require the leading lepton pT to be >= 25 GeV
  bool ZHOrthogonality = twoOSleptons && twoSFleptons && LeadpTleptons;                                   // Concatenating cuts

  // Apply both orthogonality cuts
  bool orthogonality = (GGFOrthogonality) && (!ZHOrthogonality);
  if (not Manager()->ApplyCut(orthogonality, "orthogonality"))
    return true;

  // Exactly one tight lepton
  bool oneTightLep = (tight_leptons.size() == 1);
  if (not Manager()->ApplyCut(oneTightLep, "oneTightLep"))
    return true;

  // Do Ak4 clustering
  auto Ak4result = getAk4jets(eflowTracks, eflowPhotons, eflowNeutralHadrons, tight_leptons, AK4LEP_DR);
  std::vector<fastjet::PseudoJet> Ak4jets = applyEnergyScale(Ak4result.first);
  auto ak4Pair = filter_Ak4jetsAndConstituents(Ak4jets, Ak4result.second, AK4_PT_MIN, AK4_ETA_MAX, AK4LEP_DR, tight_leptons);
  Ak4jets = ak4Pair.first;
  std::vector<std::vector<fastjet::PseudoJet>> Ak4jetConstituents = ak4Pair.second;

  //////////////////////////////////////////////
  //  Building the W and Clustering the Ak15  //
  //////////////////////////////////////////////

  // W Boson
  float const WMASS_LOW = 30;
  float const WMASS_HIGH = 130;
  float const WPT_MIN = 60;

  // Track ID
  float const TRACK_PT_MIN = 1;
  float const TRACK_ETA_MAX = 2.5;
  float const TRACK_D0_MAX = 0.05;
  float const TRACK_DZ_MAX = 0.05;
  float const TRACK_DR_MAX = 0.4;

  // Ak15 pT
  float const AK15JET_PT_MIN = 60;
  float const AK15JET_NTRACKS_MIN = 10;

  // Misc
  float const MIN_DPHI = 1.5;
  float const W_SUEP_PT_RATIO = 3.0;

  // W reconstruction
  ParticleBaseFormat recoW;
  recoW += tight_leptons.at(0).momentum();
  recoW += smearedMET;

  // Compute transverse mass: mT = sqrt(2 * pT_lep * MET * (1 - cos(dphi)))
  double lepPt = tight_leptons.at(0).pt();
  double metPt = smearedMET.pt();
  double metPhi = smearedMET.phi();
  double dphiW = computeDeltaPhi(tight_leptons.at(0).phi(), metPhi);
  double wTransverseMass = std::sqrt(2.0 * lepPt * metPt * (1.0 - std::cos(dphiW)));

  // Do Ak15 clustering
  auto Ak15result = getAk15Jets(centralTracks, tight_leptons, TRACK_PT_MIN, TRACK_ETA_MAX, TRACK_D0_MAX, TRACK_DZ_MAX, TRACK_DR_MAX);
  std::vector<fastjet::PseudoJet> Ak15Jets = Ak15result.first;
  std::vector<std::vector<fastjet::PseudoJet>> Ak15jetConstituents = Ak15result.second;

  // Do Ak15 clustering (TRACKUP)
  auto Ak15result_TRACKUP = getAk15Jets(TRACKUP, tight_leptons, TRACK_PT_MIN, TRACK_ETA_MAX, TRACK_D0_MAX, TRACK_DZ_MAX, TRACK_DR_MAX);
  std::vector<fastjet::PseudoJet> Ak15Jets_TRACKUP = Ak15result_TRACKUP.first;
  std::vector<std::vector<fastjet::PseudoJet>> Ak15jetConstituents_TRACKUP = Ak15result_TRACKUP.second;

  //////////////////////////////////////////////
  //  Second Round of Event level selections  //
  //////////////////////////////////////////////

  // Require at least one Ak15 with pT > 60 GeV
  bool oneAk15 = (Ak15Jets.size() > 0);

  bool oneAk15_TRACKUP = (Ak15Jets_TRACKUP.size() > 0);

  Manager()->ApplyCut(oneAk15, "Ak15");
  Manager()->ApplyCut(oneAk15_TRACKUP, "Ak15_TRACKUP");

  // Require at least one ak4
  bool oneAk4 = (Ak4jets.size() > 0);
  bool dROverlap = (!Ak15Jets.empty() && !Ak4jets.empty() && calculateDeltaR(Ak4jets.at(0), Ak15Jets.at(0)) < 1.5);
  if (not Manager()->ApplyCut(oneAk4 && dROverlap, "oneAk4"))
    return true;

  //////////////////////////////////////////////
  //           Filling Histograms             //
  //////////////////////////////////////////////

  // Calculating Lab and Boosted Sphericity
  double labSphericity_CENTRAL = 0.0;
  double boostedSphericity_CENTRAL = 0.0;
  double labSphericity_TRACKUP = 0.0;
  double boostedSphericity_TRACKUP = 0.0;

  if (!Ak15Jets.empty()) {
    labSphericity_CENTRAL = sphericity(Ak15jetConstituents.at(0), 1.0);
    if (Ak15jetConstituents.at(0).size() > 1) {
      auto boosted_CENTRAL = boostToSUEP(Ak15jetConstituents.at(0), Ak15Jets.at(0));
      boostedSphericity_CENTRAL = sphericity(boosted_CENTRAL, 1.0);
    } 
    else {
      boostedSphericity_CENTRAL = 0.0;
    }
  }

  if (!Ak15Jets_TRACKUP.empty()) {
    labSphericity_TRACKUP = sphericity(Ak15jetConstituents_TRACKUP.at(0), 1.0);
    if (Ak15jetConstituents_TRACKUP.at(0).size() > 1) {
      auto boosted_TRACKUP = boostToSUEP(Ak15jetConstituents_TRACKUP.at(0), Ak15Jets_TRACKUP.at(0));
      boostedSphericity_TRACKUP = sphericity(boosted_TRACKUP, 1.0);
    } 
    else {
      boostedSphericity_TRACKUP = 0.0;
    }
  }


  //////////////////////////////////////////////
  //           Writing out the CSV            //
  //////////////////////////////////////////////

  std::ofstream outFile("kinematics_output.csv", std::ios::app);

  // Write the header one time

  if (outFile.tellp() == 0) { 
    outFile << "EventNumber,lepPt,lepEta,ak15NTracks,wTransverseMass,wPt,ak4Pt,boostedSphericity_CENTRAL,metPt,looseLep,tightLep,NJets,ak15Pt,wPhi,W_SUEP_ratio\n"; 
  }

  // Event information to be stored
  // Event number
  static int eventCounter = 0;
  eventCounter++;

  outFile
      << eventCounter << ","
      << tight_leptons.at(0).pt() << ","
      << tight_leptons.at(0).eta() << ","
      << Ak15jetConstituents.at(0).size() << ","                  // size_t or unassigned long
      << wTransverseMass << ","              // double
      << recoW.pt() << ","                          // assumed float
      << Ak4jets.at(0).pt() << ","                        // assumed float
      << boostedSphericity_CENTRAL << ","            // double
      << smearedMET.pt() << ","                        // double
      << loose_leptons.size() << ","                     // int
      << tight_leptons.size() << ","                     // int
      << Ak4jets.size() << ","                        // int
      << Ak15Jets.at(0).pt() << ","                       // float
      << recoW.phi() << "\n";                       // double

  outFile.close();

  // W Histograms
  FillHistoAllRegions("wTransverseMass", wTransverseMass, regionsForHistos);
  FillHistoAllRegions("wPt", recoW.pt(), regionsForHistos);
  FillHistoAllRegions("wPhi", recoW.phi(), regionsForHistos);

  // Lepton Histograms
  FillHistoAllRegions("lepPt", tight_leptons.at(0).pt(), regionsForHistos);
  FillHistoAllRegions("lepEta", tight_leptons.at(0).eta(), regionsForHistos);
  FillHistoAllRegions("lepPhi", tight_leptons.at(0).phi(), regionsForHistos);
  FillHistoAllRegions("looseLep", loose_leptons.size(), regionsForHistos);
  FillHistoAllRegions("tightLep", tight_leptons.size(), regionsForHistos);

  // Muon or Electron Histograms
  if (tight_muons.size() > 0)
  {
    FillHistoAllRegions("muPt", tight_muons.at(0).pt(), regionsForHistos);
    FillHistoAllRegions("muEta", tight_muons.at(0).eta(), regionsForHistos);
    FillHistoAllRegions("muPhi", tight_muons.at(0).phi(), regionsForHistos);
    FillHistoAllRegions("looseMu", loose_muons.size(), regionsForHistos);
    FillHistoAllRegions("tightMu", tight_muons.size(), regionsForHistos);
  }
  else
  {
    FillHistoAllRegions("elePt", tight_electrons.at(0).pt(), regionsForHistos);
    FillHistoAllRegions("eleEta", tight_electrons.at(0).eta(), regionsForHistos);
    FillHistoAllRegions("elePhi", tight_electrons.at(0).phi(), regionsForHistos);
    FillHistoAllRegions("looseEle", loose_electrons.size(), regionsForHistos);
    FillHistoAllRegions("tightEle", tight_electrons.size(), regionsForHistos);
  }

  // Ak4 Histograms
  FillHistoAllRegions("NJets", Ak4jets.size(), regionsForHistos);
  FillHistoAllRegions("ak4Pt", Ak4jets.at(0).pt(), regionsForHistos);
  FillHistoAllRegions("ak4Eta", Ak4jets.at(0).eta(), regionsForHistos);
  FillHistoAllRegions("ak4Phi", normalizePhi(Ak4jets.at(0).phi()), regionsForHistos);
  FillHistoAllRegions("ak4NTracks", Ak4jetConstituents.at(0).size(), regionsForHistos);

  // MET Histograms
  FillHistoAllRegions("metPt", smearedMET.pt(), regionsForHistos);
  FillHistoAllRegions("metPhi", smearedMET.phi(), regionsForHistos);
  FillHistoAllRegions("dPhi_MET_Lep", computeDeltaPhi(smearedMET.phi(),tight_leptons.at(0).phi()), regionsForHistos);

  // Ak15 Histograms
  if(!Ak15Jets.empty()){
    FillHistoAllRegions("ak15Pt", Ak15Jets.at(0).pt(), regionsForHistos_CENTRAL);
    FillHistoAllRegions("ak15NTracks", Ak15jetConstituents.at(0).size(), regionsForHistos_CENTRAL);
    FillHistoAllRegions("ak15Eta", Ak15Jets.at(0).eta(), regionsForHistos_CENTRAL);
    FillHistoAllRegions("ak15Phi", normalizePhi(Ak15Jets.at(0).phi()), regionsForHistos_CENTRAL);
    FillHistoAllRegions("ak15Mass", Ak15Jets.at(0).m(), regionsForHistos_CENTRAL);
  
    FillHistoAllRegions("dPhi_MET_SUEP", computeDeltaPhi(smearedMET.phi(),normalizePhi(Ak15Jets.at(0).phi())), regionsForHistos_CENTRAL);
  
    FillHistoAllRegions("labSphericity", labSphericity_CENTRAL, regionsForHistos_CENTRAL);
    FillHistoAllRegions("boostedSphericity", boostedSphericity_CENTRAL, regionsForHistos_CENTRAL);

    if ((10 <= (Ak15jetConstituents.at(0).size()) && (Ak15jetConstituents.at(0).size()) < 20) && (0.3 <= boostedSphericity_CENTRAL && boostedSphericity_CENTRAL < 0.4))
    {
      FillHistoAllRegions("ABCD_A", Ak15jetConstituents.at(0).size(), regionsForHistos_CENTRAL);
    }
    if ((20 <= (Ak15jetConstituents.at(0).size()) && (Ak15jetConstituents.at(0).size()) < 30) && (0.3 <= boostedSphericity_CENTRAL && boostedSphericity_CENTRAL < 0.4))
    {
      FillHistoAllRegions("ABCD_B", Ak15jetConstituents.at(0).size(), regionsForHistos_CENTRAL);
    }
    if ((30 <= (Ak15jetConstituents.at(0).size()) && (Ak15jetConstituents.at(0).size()) < 200) && (0.3 <= boostedSphericity_CENTRAL && boostedSphericity_CENTRAL < 0.4))
    {
      FillHistoAllRegions("ABCD_C", Ak15jetConstituents.at(0).size(), regionsForHistos_CENTRAL);
    }

    if ((10 <= (Ak15jetConstituents.at(0).size()) && (Ak15jetConstituents.at(0).size()) < 20) && (0.4 <= boostedSphericity_CENTRAL && boostedSphericity_CENTRAL < 0.5))
    {
      FillHistoAllRegions("ABCD_D", Ak15jetConstituents.at(0).size(), regionsForHistos_CENTRAL);
    }
    if ((20 <= (Ak15jetConstituents.at(0).size()) && (Ak15jetConstituents.at(0).size()) < 30) && (0.4 <= boostedSphericity_CENTRAL && boostedSphericity_CENTRAL < 0.5))
    {
      FillHistoAllRegions("ABCD_E", Ak15jetConstituents.at(0).size(), regionsForHistos_CENTRAL);
    }
    if ((30 <= (Ak15jetConstituents.at(0).size()) && (Ak15jetConstituents.at(0).size()) < 40) && (0.4 <= boostedSphericity_CENTRAL && boostedSphericity_CENTRAL < 0.5))
    {
      FillHistoAllRegions("ABCD_F0", Ak15jetConstituents.at(0).size(), regionsForHistos_CENTRAL);
    }
    if ((40 <= (Ak15jetConstituents.at(0).size()) && (Ak15jetConstituents.at(0).size()) < 50) && (0.4 <= boostedSphericity_CENTRAL && boostedSphericity_CENTRAL < 0.5))
    {
      FillHistoAllRegions("ABCD_F1", Ak15jetConstituents.at(0).size(), regionsForHistos_CENTRAL);
    }
    if ((50 <= (Ak15jetConstituents.at(0).size()) && (Ak15jetConstituents.at(0).size()) < 60) && (0.4 <= boostedSphericity_CENTRAL && boostedSphericity_CENTRAL < 0.5))
    {
      FillHistoAllRegions("ABCD_F2", Ak15jetConstituents.at(0).size(), regionsForHistos_CENTRAL);
    }
    if ((60 <= (Ak15jetConstituents.at(0).size()) && (Ak15jetConstituents.at(0).size()) < 80) && (0.4 <= boostedSphericity_CENTRAL && boostedSphericity_CENTRAL < 0.5))
    {
      FillHistoAllRegions("ABCD_F3", Ak15jetConstituents.at(0).size(), regionsForHistos_CENTRAL);
    }
    if ((80 <= (Ak15jetConstituents.at(0).size()) && (Ak15jetConstituents.at(0).size()) < 200) && (0.4 <= boostedSphericity_CENTRAL && boostedSphericity_CENTRAL < 0.5))
    {
      FillHistoAllRegions("ABCD_F4", Ak15jetConstituents.at(0).size(), regionsForHistos_CENTRAL);
    }

    if ((10 <= (Ak15jetConstituents.at(0).size()) && (Ak15jetConstituents.at(0).size()) < 20) && (0.5 <= boostedSphericity_CENTRAL && boostedSphericity_CENTRAL <= 1.0))
    {
      FillHistoAllRegions("ABCD_G", Ak15jetConstituents.at(0).size(), regionsForHistos_CENTRAL);
    }
    if ((20 <= (Ak15jetConstituents.at(0).size()) && (Ak15jetConstituents.at(0).size()) < 30) && (0.5 <= boostedSphericity_CENTRAL && boostedSphericity_CENTRAL <= 1.0))
    {
      FillHistoAllRegions("ABCD_H", Ak15jetConstituents.at(0).size(), regionsForHistos_CENTRAL);
    }
    if ((30 <= (Ak15jetConstituents.at(0).size()) && (Ak15jetConstituents.at(0).size()) < 40) && (0.5 <= boostedSphericity_CENTRAL && boostedSphericity_CENTRAL <= 1.0))
    {
      FillHistoAllRegions("ABCD_SR0", Ak15jetConstituents.at(0).size(), regionsForHistos_CENTRAL);
    }
    if ((40 <= (Ak15jetConstituents.at(0).size()) && (Ak15jetConstituents.at(0).size()) < 50) && (0.5 <= boostedSphericity_CENTRAL && boostedSphericity_CENTRAL <= 1.0))
    {
      FillHistoAllRegions("ABCD_SR1", Ak15jetConstituents.at(0).size(), regionsForHistos_CENTRAL);
    }
    if ((50 <= (Ak15jetConstituents.at(0).size()) && (Ak15jetConstituents.at(0).size()) < 60) && (0.5 <= boostedSphericity_CENTRAL && boostedSphericity_CENTRAL <= 1.0))
    {
      FillHistoAllRegions("ABCD_SR2", Ak15jetConstituents.at(0).size(), regionsForHistos_CENTRAL);
    }
    if ((60 <= (Ak15jetConstituents.at(0).size()) && (Ak15jetConstituents.at(0).size()) < 80) && (0.5 <= boostedSphericity_CENTRAL && boostedSphericity_CENTRAL <= 1.0))
    {
      FillHistoAllRegions("ABCD_SR3", Ak15jetConstituents.at(0).size(), regionsForHistos_CENTRAL);
    }
    if ((80 <= (Ak15jetConstituents.at(0).size()) && (Ak15jetConstituents.at(0).size()) < 200) && (0.5 <= boostedSphericity_CENTRAL && boostedSphericity_CENTRAL <= 1.0))
    {
      FillHistoAllRegions("ABCD_SR4", Ak15jetConstituents.at(0).size(), regionsForHistos_CENTRAL);
    }

    // ZH binnning
    if ((Ak15jetConstituents.at(0).size() >= 21.0) && (Ak4jets.at(0).pt() <= 135.0))
    {
      FillHistoAllRegions("ABCD_A_ZH", Ak15jetConstituents.at(0).size(), regionsForHistos_CENTRAL);
    }
    if ((Ak15jetConstituents.at(0).size() >= 14.0) && (Ak15jetConstituents.at(0).size() < 21.0) && (Ak4jets.at(0).pt() <= 135.0))
    {
      FillHistoAllRegions("ABCD_C1_ZH", Ak15jetConstituents.at(0).size(), regionsForHistos_CENTRAL);
    }
    if ((Ak15jetConstituents.at(0).size() < 14.0) && (Ak4jets.at(0).pt() <= 135.0))
    {
      FillHistoAllRegions("ABCD_C2_ZH", Ak15jetConstituents.at(0).size(), regionsForHistos_CENTRAL);
    }

    if ((Ak15jetConstituents.at(0).size() >= 21.0) && (Ak4jets.at(0).pt() <= 220.0) && (Ak4jets.at(0).pt() > 135.0))
    {
      FillHistoAllRegions("ABCD_B1_ZH", Ak15jetConstituents.at(0).size(), regionsForHistos_CENTRAL);
    }
    if ((Ak15jetConstituents.at(0).size() >= 14.0) && (Ak15jetConstituents.at(0).size() < 21.0) && (Ak4jets.at(0).pt() <= 220.0) && (Ak4jets.at(0).pt() > 135.0))
    {
      FillHistoAllRegions("ABCD_D1_ZH", Ak15jetConstituents.at(0).size(), regionsForHistos_CENTRAL);
    }
    if ((Ak15jetConstituents.at(0).size() < 14.0) && (Ak4jets.at(0).pt() <= 220.0) && (Ak4jets.at(0).pt() > 135.0))
    {
      FillHistoAllRegions("ABCD_D2_ZH", Ak15jetConstituents.at(0).size(), regionsForHistos_CENTRAL);
    }

    if ((Ak15jetConstituents.at(0).size() >= 21.0) && (Ak4jets.at(0).pt() > 220.0))
    {
      FillHistoAllRegions("ABCD_B2_ZH", Ak15jetConstituents.at(0).size(), regionsForHistos_CENTRAL);
    }
    if ((Ak15jetConstituents.at(0).size() >= 14.0) && (Ak15jetConstituents.at(0).size() < 21.0) && (Ak4jets.at(0).pt() > 220.0))
    {
      FillHistoAllRegions("ABCD_E1_ZH", Ak15jetConstituents.at(0).size(), regionsForHistos_CENTRAL);
    }
    if ((Ak15jetConstituents.at(0).size() < 14.0) && (Ak4jets.at(0).pt() > 220.0))
    {
      FillHistoAllRegions("ABCD_E2_ZH", Ak15jetConstituents.at(0).size(), regionsForHistos_CENTRAL);
    }
  }

  // Ak15 Histograms (TRACKUP)
  if(!Ak15Jets_TRACKUP.empty()){
    FillHistoAllRegions("ak15Pt", Ak15Jets_TRACKUP.at(0).pt(), regionsForHistos_TRACKUP);
    FillHistoAllRegions("ak15NTracks", Ak15jetConstituents_TRACKUP.at(0).size(), regionsForHistos_TRACKUP);
    FillHistoAllRegions("ak15Eta", Ak15Jets_TRACKUP.at(0).eta(), regionsForHistos_TRACKUP);
    FillHistoAllRegions("ak15Phi", normalizePhi(Ak15Jets_TRACKUP.at(0).phi()), regionsForHistos_TRACKUP);
    FillHistoAllRegions("ak15Mass", Ak15Jets_TRACKUP.at(0).m(), regionsForHistos_TRACKUP);
  
    FillHistoAllRegions("dPhi_MET_SUEP", computeDeltaPhi(smearedMET.phi(),normalizePhi(Ak15Jets_TRACKUP.at(0).phi())), regionsForHistos_TRACKUP);
  
    FillHistoAllRegions("labSphericity", labSphericity_TRACKUP, regionsForHistos_TRACKUP);
    FillHistoAllRegions("boostedSphericity", boostedSphericity_TRACKUP, regionsForHistos_TRACKUP);

    if ((10 <= (Ak15jetConstituents_TRACKUP.at(0).size()) && (Ak15jetConstituents_TRACKUP.at(0).size()) < 20) && (0.3 <= boostedSphericity_TRACKUP && boostedSphericity_TRACKUP < 0.4))
    {
      FillHistoAllRegions("ABCD_A", Ak15jetConstituents_TRACKUP.at(0).size(), regionsForHistos_TRACKUP);
    }
    if ((20 <= (Ak15jetConstituents_TRACKUP.at(0).size()) && (Ak15jetConstituents_TRACKUP.at(0).size()) < 30) && (0.3 <= boostedSphericity_TRACKUP && boostedSphericity_TRACKUP < 0.4))
    {
      FillHistoAllRegions("ABCD_B", Ak15jetConstituents_TRACKUP.at(0).size(), regionsForHistos_TRACKUP);
    }
    if ((30 <= (Ak15jetConstituents_TRACKUP.at(0).size()) && (Ak15jetConstituents_TRACKUP.at(0).size()) < 200) && (0.3 <= boostedSphericity_TRACKUP && boostedSphericity_TRACKUP < 0.4))
    {
      FillHistoAllRegions("ABCD_C", Ak15jetConstituents_TRACKUP.at(0).size(), regionsForHistos_TRACKUP);
    }

    if ((10 <= (Ak15jetConstituents_TRACKUP.at(0).size()) && (Ak15jetConstituents_TRACKUP.at(0).size()) < 20) && (0.4 <= boostedSphericity_TRACKUP && boostedSphericity_TRACKUP < 0.5))
    {
      FillHistoAllRegions("ABCD_D", Ak15jetConstituents_TRACKUP.at(0).size(), regionsForHistos_TRACKUP);
    }
    if ((20 <= (Ak15jetConstituents_TRACKUP.at(0).size()) && (Ak15jetConstituents_TRACKUP.at(0).size()) < 30) && (0.4 <= boostedSphericity_TRACKUP && boostedSphericity_TRACKUP < 0.5))
    {
      FillHistoAllRegions("ABCD_E", Ak15jetConstituents_TRACKUP.at(0).size(), regionsForHistos_TRACKUP);
    }
    if ((30 <= (Ak15jetConstituents_TRACKUP.at(0).size()) && (Ak15jetConstituents_TRACKUP.at(0).size()) < 40) && (0.4 <= boostedSphericity_TRACKUP && boostedSphericity_TRACKUP < 0.5))
    {
      FillHistoAllRegions("ABCD_F0", Ak15jetConstituents_TRACKUP.at(0).size(), regionsForHistos_TRACKUP);
    }
    if ((40 <= (Ak15jetConstituents_TRACKUP.at(0).size()) && (Ak15jetConstituents_TRACKUP.at(0).size()) < 50) && (0.4 <= boostedSphericity_TRACKUP && boostedSphericity_TRACKUP < 0.5))
    {
      FillHistoAllRegions("ABCD_F1", Ak15jetConstituents_TRACKUP.at(0).size(), regionsForHistos_TRACKUP);
    }
    if ((50 <= (Ak15jetConstituents_TRACKUP.at(0).size()) && (Ak15jetConstituents_TRACKUP.at(0).size()) < 60) && (0.4 <= boostedSphericity_TRACKUP && boostedSphericity_TRACKUP < 0.5))
    {
      FillHistoAllRegions("ABCD_F2", Ak15jetConstituents_TRACKUP.at(0).size(), regionsForHistos_TRACKUP);
    }
    if ((60 <= (Ak15jetConstituents_TRACKUP.at(0).size()) && (Ak15jetConstituents_TRACKUP.at(0).size()) < 80) && (0.4 <= boostedSphericity_TRACKUP && boostedSphericity_TRACKUP < 0.5))
    {
      FillHistoAllRegions("ABCD_F3", Ak15jetConstituents_TRACKUP.at(0).size(), regionsForHistos_TRACKUP);
    }
    if ((80 <= (Ak15jetConstituents_TRACKUP.at(0).size()) && (Ak15jetConstituents_TRACKUP.at(0).size()) < 200) && (0.4 <= boostedSphericity_TRACKUP && boostedSphericity_TRACKUP < 0.5))
    {
      FillHistoAllRegions("ABCD_F4", Ak15jetConstituents_TRACKUP.at(0).size(), regionsForHistos_TRACKUP);
    }

    if ((10 <= (Ak15jetConstituents_TRACKUP.at(0).size()) && (Ak15jetConstituents_TRACKUP.at(0).size()) < 20) && (0.5 <= boostedSphericity_TRACKUP && boostedSphericity_TRACKUP <= 1.0))
    {
      FillHistoAllRegions("ABCD_G", Ak15jetConstituents_TRACKUP.at(0).size(), regionsForHistos_TRACKUP);
    }
    if ((20 <= (Ak15jetConstituents_TRACKUP.at(0).size()) && (Ak15jetConstituents_TRACKUP.at(0).size()) < 30) && (0.5 <= boostedSphericity_TRACKUP && boostedSphericity_TRACKUP <= 1.0))
    {
      FillHistoAllRegions("ABCD_H", Ak15jetConstituents_TRACKUP.at(0).size(), regionsForHistos_TRACKUP);
    }
    if ((30 <= (Ak15jetConstituents_TRACKUP.at(0).size()) && (Ak15jetConstituents_TRACKUP.at(0).size()) < 40) && (0.5 <= boostedSphericity_TRACKUP && boostedSphericity_TRACKUP <= 1.0))
    {
      FillHistoAllRegions("ABCD_SR0", Ak15jetConstituents_TRACKUP.at(0).size(), regionsForHistos_TRACKUP);
    }
    if ((40 <= (Ak15jetConstituents_TRACKUP.at(0).size()) && (Ak15jetConstituents_TRACKUP.at(0).size()) < 50) && (0.5 <= boostedSphericity_TRACKUP && boostedSphericity_TRACKUP <= 1.0))
    {
      FillHistoAllRegions("ABCD_SR1", Ak15jetConstituents_TRACKUP.at(0).size(), regionsForHistos_TRACKUP);
    }
    if ((50 <= (Ak15jetConstituents_TRACKUP.at(0).size()) && (Ak15jetConstituents_TRACKUP.at(0).size()) < 60) && (0.5 <= boostedSphericity_TRACKUP && boostedSphericity_TRACKUP <= 1.0))
    {
      FillHistoAllRegions("ABCD_SR2", Ak15jetConstituents_TRACKUP.at(0).size(), regionsForHistos_TRACKUP);
    }
    if ((60 <= (Ak15jetConstituents_TRACKUP.at(0).size()) && (Ak15jetConstituents_TRACKUP.at(0).size()) < 80) && (0.5 <= boostedSphericity_TRACKUP && boostedSphericity_TRACKUP <= 1.0))
    {
      FillHistoAllRegions("ABCD_SR3", Ak15jetConstituents_TRACKUP.at(0).size(), regionsForHistos_TRACKUP);
    }
    if ((80 <= (Ak15jetConstituents_TRACKUP.at(0).size()) && (Ak15jetConstituents_TRACKUP.at(0).size()) < 200) && (0.5 <= boostedSphericity_TRACKUP && boostedSphericity_TRACKUP <= 1.0))
    {
      FillHistoAllRegions("ABCD_SR4", Ak15jetConstituents_TRACKUP.at(0).size(), regionsForHistos_TRACKUP);
    }

  }



  return true;
}

///////////////////////////////////////////////////////////////
//                        Finalize                           //
//    Function called one time at the end of the analysis    //
///////////////////////////////////////////////////////////////

void user::Finalize(const SampleFormat &summary, const std::vector<SampleFormat> &files) {}

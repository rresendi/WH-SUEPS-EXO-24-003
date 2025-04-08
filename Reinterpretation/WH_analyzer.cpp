#include "SampleAnalyzer/User/Analyzer/user.h"
using namespace MA5;
using namespace std;
#include <Eigen/Dense>

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

bool iso(const RecLeptonFormat &lepton,
         const std::vector<RecTrackFormat> &tracks,
         const std::vector<RecParticleFormat> &eflowPhotons,
         const std::vector<RecParticleFormat> &eflowNeutralHadrons,
         double iso_minpt,
         double deltaRmax,
         std::string iso = "")
{

    double ratiomax = 0.0;

    // to be played with
    if (iso == "WP90")
    {
        ratiomax = 0.12;
    }
    else if (iso == "WP80")
    {
        ratiomax = 0.10;
    }
    else if (iso == "pfIso2")
    {
        ratiomax = 0.25;
    }
    else if (iso == "pfIso5")
    {
        ratiomax = 0.2;
    }

    double lep_pt = lepton.pt();
    double totalpt = 0.0;

    for (const auto &track : tracks)
    {
        double dr = lepton.dr(track);
        double pt = track.pt();
        // int pdgid = std::abs(track.pdgid()); // Bug in MadAnalysis before we can use this
        if (dr < 0.001)
        {
            continue;
        }

        // if (pdgid == 11 || pdgid == 13)
        //{
        //     continue;
        // }

        if (dr < deltaRmax && pt > iso_minpt)
        {
            totalpt += pt;
        }
    }

    for (const auto &photon : eflowPhotons)
    {
        double dr = lepton.dr(photon);
        double pt = photon.pt();
        if (dr < deltaRmax && pt > iso_minpt)
        {
            totalpt += pt;
        }
    }

    for (const auto &neutral : eflowNeutralHadrons)
    {
        double dr = lepton.dr(neutral);
        double pt = neutral.pt();
        if (dr < deltaRmax && pt > iso_minpt)
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

        // Loose Object Selections
        if (selection == "loose")
        {

            if (fabs(obj.eta()) > etamax)
                continue;

            if (fabs(obj.d0()) > d0)
                continue;
        }

        // Loose or tight
        if (obj.pt() < ptmin)
            continue;

        if (fabs(obj.dz()) > dz)
            continue;

        // Iso selection: use "pfIso2" for loose and "pfIso5" for tight
        string iso_option = (selection == "loose") ? "pfIso2" : "pfIso5";
        if (!iso(obj, tracks, eflowPhotons, eflowNeutralHadrons, iso_pTMin, iso_dRMax, iso_option))
            continue;

        filtered.push_back(obj);
    }

    return filtered;
}

vector<RecLeptonFormat> filter_electrons(const vector<RecLeptonFormat> &objects,
                                         float ptmin,
                                         float etamax,
                                         float iso_pTMin,
                                         float iso_dRMax,
                                         const vector<RecTrackFormat> &tracks,
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

        // Loose or tight
        if (obj.pt() < ptmin)
            continue;

        // Loose Object Selections
        if (selection == "loose")
        {
            if (fabs(obj.eta()) > etamax)
                continue;

            // Exclude crack region: electrons with 1.444 < |eta| < 1.566
            if ((fabs(obj.eta()) > 1.444) && (fabs(obj.eta()) < 1.566))
                continue;

            if (fabs(obj.d0()) > (0.05 + 0.05 * (fabs(obj.eta()) > 1.479)))
                continue;

            if (fabs(obj.dz()) > (0.1 + 0.1 * (fabs(obj.eta()) > 1.479)))
                continue;
        }

        string iso_option = (selection == "loose") ? "WP90" : "WP80";
        if (!iso(obj, tracks, eflowPhotons, eflowNeutralHadrons, iso_pTMin, iso_dRMax, iso_option))
        {
            continue;
        }

        filtered.push_back(obj);
    }
    return filtered;
}

vector<RecJetFormat> filter_Ak4jets(vector<RecJetFormat> jets,
                                    float ptmin,
                                    float etamax,
                                    float deltaRmin,
                                    const vector<RecLeptonFormat> &leptons)
{
    // Helper function to select jets passing pt, eta, and ΔR selections
    vector<RecJetFormat> filtered;

    for (auto &jet : jets)
    {

        // Object Selections
        if (jet.pt() < ptmin)
        {
            continue;
        }
        if (fabs(jet.eta()) > etamax)
        {
            continue;
        }

        // ΔR selection (jet-lepton separation)
        bool passesDeltaRCut = true;

        // Loop over each lepton to check if any are too close to the jet
        for (const auto &lepton : leptons)
        {
            // Calculate ΔR between the jet and the lepton
            if (jet.dr(lepton) <= deltaRmin)
            {
                passesDeltaRCut = false;
                break; // If one lepton is too close, reject this jet
            }
        }

        // If the jet passes the ΔR cut, add it to the filtered list
        if (passesDeltaRCut)
        {
            filtered.push_back(jet);
        }
    }

    return filtered;
}

bool compareBypT(const RecLeptonFormat &a, const RecLeptonFormat &b)
{
    return a.pt() > b.pt(); // Sort by pT in descending order
}

bool compareBypTJets(const RecJetFormat &a, const RecJetFormat &b)
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

// Function to sort jet collections by pT
void sortJetCollections(std::vector<RecJetFormat> &Ak4jets)
{
    if (!Ak4jets.empty())
    {
        std::sort(Ak4jets.begin(), Ak4jets.end(), compareBypTJets);
    }
}

// Veto b-tagged jets

bool BTagVeto(vector<RecJetFormat> jets)
{
    // Returns true if a b jet is found
    for (auto jet : jets)
    {
        if (jet.btag())
        {
            return false;
        }
    }
    return true;
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
            track.dr(leptons.at(0)) >= dr_cut)

        // Creating the particles (constituents) which makeup our jet
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

// calculating deltaR between the reconstructed jet and the simulated jet

double calculateDeltaR(const RecJetFormat &rec_jet, const fastjet::PseudoJet &pseudo_jet)
{
    double rec_eta = rec_jet.eta();
    double rec_phi = rec_jet.phi();

    double pseudo_eta = pseudo_jet.eta();
    double pseudo_phi = pseudo_jet.phi();

    // Map pseudo_phi from [0, 2pi] to [-pi, pi] for consistency
    if (pseudo_phi > M_PI)
    {
        pseudo_phi -= 2 * M_PI;
    }

    // Calculate delta eta and delta phi
    double delta_eta = rec_eta - pseudo_eta;
    double delta_phi = std::atan2(std::sin(rec_phi - pseudo_phi), std::cos(rec_phi - pseudo_phi));

    return std::sqrt(delta_eta * delta_eta + delta_phi * delta_phi);
}

bool dPhiCut(float phi1, float phi2, float MIN_DPHI)
{
    float dphi = phi1 - phi2;
    if (dphi > M_PI)
    {
        dphi -= 2 * M_PI;
    }
    else if (dphi < -M_PI)
    {
        dphi += 2 * M_PI;
    }
    return (std::fabs(dphi) > MIN_DPHI);
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
    PHYSICS->recConfig().UseDeltaRIsolation(0.5);

    // ===== Signal region ===== //
    Manager()->AddRegionSelection("SR");

    // ===== Selections ===== //
    Manager()->AddCut("orthogonality");
    Manager()->AddCut("oneTightLep");
    Manager()->AddCut("oneAk4");
    Manager()->AddCut("firstMETCut");
    Manager()->AddCut("Ak15");
    Manager()->AddCut("secondMETCut");
    Manager()->AddCut("wPtCut");
    Manager()->AddCut("onShellW");
    Manager()->AddCut("noBTag");
    Manager()->AddCut("dPhi_W_SUEP");
    Manager()->AddCut("dPhi_MET_SUEP");
    Manager()->AddCut("dPhi_Lep_SUEP");
    Manager()->AddCut("Ak4_Ak15_Overlap");
    Manager()->AddCut("W_SUEP_ratio");
    Manager()->AddCut("dPhi_MET_Ak4");

    // ===== Histograms ===== //

    // W Histograms
    Manager()->AddHisto("wMass", 60, 0.0, 140.0);
    Manager()->AddHisto("wPt", 300, 0, 300);
    Manager()->AddHisto("wEta", 40, -3.14, 3.14);
    Manager()->AddHisto("wPhi", 40, -3.14, 3.14);

    // Lepton Histograms
    Manager()->AddHisto("lepPt", 300, 0, 300);
    Manager()->AddHisto("lepEta", 40, -3.14, 3.14);
    Manager()->AddHisto("lepPhi", 40, -3.14, 3.14);
    Manager()->AddHisto("looseLep", 10, 0, 10);
    Manager()->AddHisto("tightLep", 10, 0, 10);

    // Muon Histograms
    Manager()->AddHisto("muPt", 300, 0, 300);
    Manager()->AddHisto("muEta", 40, -3.14, 3.14);
    Manager()->AddHisto("muPhi", 40, -3.14, 3.14);
    Manager()->AddHisto("looseMu", 10, 0, 10);
    Manager()->AddHisto("tightMu", 10, 0, 10);

    // Electron Histograms
    Manager()->AddHisto("elePt", 300, 0, 300);
    Manager()->AddHisto("eleEta", 40, -3.14, 3.14);
    Manager()->AddHisto("elePhi", 40, -3.14, 3.14);
    Manager()->AddHisto("looseEle", 10, 0, 10);
    Manager()->AddHisto("tightEle", 10, 0, 10);

    // Ak4 Histograms
    Manager()->AddHisto("NJets", 10, 0.0, 10.0);
    Manager()->AddHisto("ak4Pt", 300, 0, 300);
    Manager()->AddHisto("ak4Eta", 40, -3.14, 3.14);
    Manager()->AddHisto("ak4Phi", 40, -3.14, 3.14);
    Manager()->AddHisto("ak4NTracks", 100, 0.0, 100.0);

    // Ak15 Histograms
    Manager()->AddHisto("ak15Pt", 300, 0, 300);
    Manager()->AddHisto("ak15Eta", 40, -3.14, 3.14);
    Manager()->AddHisto("ak15Phi", 40, -3.14, 3.14);
    Manager()->AddHisto("ak15NTracks", 100, 0.0, 100.0);
    Manager()->AddHisto("ak15Mass", 400, 0, 400);

    // Sphericity Histograms
    Manager()->AddHisto("labSphericity", 50, 0.0, 1.0);
    Manager()->AddHisto("boostedSphericity", 50, 0.0, 1.0);

    // MET Histograms
    Manager()->AddHisto("metPt", 300, 0, 300);
    Manager()->AddHisto("metPhi", 40, -3.14, 3.14);
    Manager()->AddHisto("metEta", 40, -3.14, 3.14);

    // ABCD Histograms
    Manager()->AddHisto("ABCD_A", 10, 10.0, 20.0);
    Manager()->AddHisto("ABCD_B", 10, 20.0, 30.0);
    Manager()->AddHisto("ABCD_C", 170, 30.0, 200.0);
    Manager()->AddHisto("ABCD_D", 10, 10.0, 20.0);
    Manager()->AddHisto("ABCD_E", 10, 20.0, 30.0);
    Manager()->AddHisto("ABCD_F0", 10, 30.0, 40.0);
    Manager()->AddHisto("ABCD_F1", 10, 40.0, 50.0);
    Manager()->AddHisto("ABCD_F2", 10, 50.0, 60.0);
    Manager()->AddHisto("ABCD_F3", 20, 60.0, 80.0);
    Manager()->AddHisto("ABCD_F4", 120, 80.0, 200.0);
    Manager()->AddHisto("ABCD_G", 10, 10.0, 20.0);
    Manager()->AddHisto("ABCD_H", 10, 20.0, 30.0);
    Manager()->AddHisto("ABCD_SR0", 10, 30.0, 40.0);
    Manager()->AddHisto("ABCD_SR1", 10, 40.0, 50.0);
    Manager()->AddHisto("ABCD_SR2", 10, 50.0, 60.0);
    Manager()->AddHisto("ABCD_SR3", 20, 60.0, 80.0);
    Manager()->AddHisto("ABCD_SR4", 20, 120.0, 200.0);

    // everything runs smoothly //
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

    // DEFINING OBJECT CUTS

    // Isolation Eflow Collections
    std::vector<RecTrackFormat> tracks = event.rec()->tracks();
    std::vector<RecParticleFormat> eflowPhotons = event.rec()->EFlowPhotons();
    std::vector<RecParticleFormat> eflowNeutralHadrons = event.rec()->EFlowNeutralHadrons();

    // ELECTRONS
    float const LOOSE_ELECTRON_PT_MIN = 15;
    float const TIGHT_ELECTRON_PT_MIN = 35;
    float const ELECTRON_ETA_MAX = 2.5;
    float const ELECTRON_ISO_PT_MIN = 0.5;
    float const ELECTRON_ISO_DR_MAX = 0.3;

    // MUONS
    float const LOOSE_MUON_PT_MIN = 10;
    float const TIGHT_MUON_PT_MIN = 30;
    float const MUON_ETA_MAX = 2.4;
    float const LOOSE_MUON_DZ = 0.1;
    float const TIGHT_MUON_DZ = 0.05;
    float const MUON_D0 = 0.02;
    float const MUON_ISO_PT_MIN = 0.1;
    float const MUON_ISO_DR_MAX = 0.4;

    // Ak4jets
    float const AK4_PT_MIN = 30;
    float const AK4_ETA_MAX = 2.4;
    float const AK4LEP_DR = 0.4;

    //////////////////////////////////////////////
    //  Applying Base Lepton Object Selections  //
    //////////////////////////////////////////////

    // Electron Collections
    vector<RecLeptonFormat> loose_electrons = filter_electrons(event.rec()->electrons(), LOOSE_ELECTRON_PT_MIN, ELECTRON_ETA_MAX, ELECTRON_ISO_PT_MIN, ELECTRON_ISO_DR_MAX, tracks, eflowPhotons, eflowNeutralHadrons, "loose");
    vector<RecLeptonFormat> loose_posElectrons = filter_electrons(loose_electrons, LOOSE_ELECTRON_PT_MIN, ELECTRON_ETA_MAX, ELECTRON_ISO_PT_MIN, ELECTRON_ISO_DR_MAX, tracks, eflowPhotons, eflowNeutralHadrons, "loose", "+");
    vector<RecLeptonFormat> loose_negElectrons = filter_electrons(loose_electrons, LOOSE_ELECTRON_PT_MIN, ELECTRON_ETA_MAX, ELECTRON_ISO_PT_MIN, ELECTRON_ISO_DR_MAX, tracks, eflowPhotons, eflowNeutralHadrons, "loose", "-");

    vector<RecLeptonFormat> tight_electrons = filter_electrons(loose_electrons, TIGHT_ELECTRON_PT_MIN, ELECTRON_ETA_MAX, ELECTRON_ISO_PT_MIN, ELECTRON_ISO_DR_MAX, tracks, eflowPhotons, eflowNeutralHadrons, "tight");
    vector<RecLeptonFormat> tight_posElectrons = filter_electrons(tight_electrons, TIGHT_ELECTRON_PT_MIN, ELECTRON_ETA_MAX, ELECTRON_ISO_PT_MIN, ELECTRON_ISO_DR_MAX, tracks, eflowPhotons, eflowNeutralHadrons, "tight", "+");
    vector<RecLeptonFormat> tight_negElectrons = filter_electrons(tight_electrons, TIGHT_ELECTRON_PT_MIN, ELECTRON_ETA_MAX, ELECTRON_ISO_PT_MIN, ELECTRON_ISO_DR_MAX, tracks, eflowPhotons, eflowNeutralHadrons, "tight", "-");

    // Muon Collections
    vector<RecLeptonFormat> loose_muons = filter_muons(event.rec()->muons(), LOOSE_MUON_PT_MIN, MUON_ETA_MAX, MUON_D0, LOOSE_MUON_DZ, MUON_ISO_PT_MIN, MUON_ISO_DR_MAX, tracks, eflowPhotons, eflowNeutralHadrons, "loose");
    vector<RecLeptonFormat> loose_posMuons = filter_muons(loose_muons, LOOSE_MUON_PT_MIN, MUON_ETA_MAX, MUON_D0, LOOSE_MUON_DZ, MUON_ISO_PT_MIN, MUON_ISO_DR_MAX, tracks, eflowPhotons, eflowNeutralHadrons, "loose", "+");
    vector<RecLeptonFormat> loose_negMuons = filter_muons(loose_muons, LOOSE_MUON_PT_MIN, MUON_ETA_MAX, MUON_D0, LOOSE_MUON_DZ, MUON_ISO_PT_MIN, MUON_ISO_DR_MAX, tracks, eflowPhotons, eflowNeutralHadrons, "loose", "-");

    vector<RecLeptonFormat> tight_muons = filter_muons(loose_muons, TIGHT_MUON_PT_MIN, MUON_ETA_MAX, MUON_D0, TIGHT_MUON_DZ, MUON_ISO_PT_MIN, MUON_ISO_DR_MAX, tracks, eflowPhotons, eflowNeutralHadrons, "tight");
    vector<RecLeptonFormat> tight_posMuons = filter_muons(tight_muons, TIGHT_MUON_PT_MIN, MUON_ETA_MAX, MUON_D0, TIGHT_MUON_DZ, MUON_ISO_PT_MIN, MUON_ISO_DR_MAX, tracks, eflowPhotons, eflowNeutralHadrons, "tight", "+");
    vector<RecLeptonFormat> tight_negMuons = filter_muons(tight_muons, TIGHT_MUON_PT_MIN, MUON_ETA_MAX, MUON_D0, TIGHT_MUON_DZ, MUON_ISO_PT_MIN, MUON_ISO_DR_MAX, tracks, eflowPhotons, eflowNeutralHadrons, "tight", "-");

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
    bool GGFOrthogonality = (loose_leptons.size() == 0);

    // Orthogonality to ZH: Remove events with a pair of loose OSSF leptons
    bool twoOSleptons = (loose_posLeptons.size() == 1 && loose_negLeptons.size() == 1);                     // Require exactly two opposite-sign leptons (either muons or electrons)
    bool twoSFleptons = (loose_muons.size() == 2 || loose_electrons.size() == 2);                           // Require exactly two muons or exactly two electrons
    bool LeadpTleptons = (loose_leptons.size() > 0) && (loose_leptons.at(0).pt() >= ORTHOG_LEAD_LEPTON_PT); // Require the leading lepton pT to be >= 25 GeV
    bool ZHOrthogonality = twoOSleptons && twoSFleptons && LeadpTleptons;                                   // Concatenating cuts

    // Apply both orthogonality cuts
    bool orthogonality = (!GGFOrthogonality) && (!ZHOrthogonality);
    if (not Manager()->ApplyCut(orthogonality, "orthogonality"))
        return true;

    // Exactly one tight lepton
    bool oneTightLep = (tight_leptons.size() == 1);
    if (not Manager()->ApplyCut(oneTightLep, "oneTightLep"))
        return true;

    // Define the Ak4 jet collection with our tight lepton for cleaning
    vector<RecJetFormat> Ak4jets = filter_Ak4jets(event.rec()->jets(), AK4_PT_MIN, AK4_ETA_MAX, AK4LEP_DR, tight_leptons);

    // Sorting Ak4 jets
    sortJetCollections(Ak4jets);

    // Require at least one ak4
    bool oneAk4 = (Ak4jets.size() > 0);
    if (not Manager()->ApplyCut(oneAk4, "oneAk4"))
        return true;

    // First MET pT Cut
    bool firstMETCut = (event.rec()->MET().pt() > FIRST_MET_PT);
    if (!Manager()->ApplyCut(firstMETCut, "firstMETCut"))
        return true;

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

    // Misc
    float const MIN_DPHI = 1.5;
    float const W_SUEP_PT_RATIO = 3.0;

    // W reconstruction
    ParticleBaseFormat recoW;
    recoW += tight_leptons.at(0).momentum();
    recoW += event.rec()->MET().momentum();

    // Do Ak15 clustering
    auto Ak15result = getAk15Jets(event.rec()->tracks(), tight_leptons, TRACK_PT_MIN, TRACK_ETA_MAX, TRACK_D0_MAX, TRACK_DZ_MAX, TRACK_DR_MAX);

    std::vector<fastjet::PseudoJet> Ak15Jets = Ak15result.first;
    std::vector<std::vector<fastjet::PseudoJet>> Ak15jetConstituents = Ak15result.second;

    //////////////////////////////////////////////
    //  Second Round of Event level selections  //
    //////////////////////////////////////////////

    // Require at least one Ak15 with pT > 60 GeV
    bool oneAk15 = (Ak15Jets.size() > 0);
    bool minAk15pT = false;
    if (oneAk15)
    {
        minAk15pT = (Ak15Jets.at(0).pt() > AK15JET_PT_MIN);
    }
    bool Ak15 = (oneAk15) && (minAk15pT);
    if (not Manager()->ApplyCut(Ak15, "Ak15"))
        return true;

    // Second MET pT Cut
    bool secondMETCut = (event.rec()->MET().pt() > SECOND_MET_PT);
    if (!Manager()->ApplyCut(secondMETCut, "secondMETCut"))
        return true;

    // Require W pT
    bool wPtCut = (recoW.pt() > WPT_MIN);
    if (not Manager()->ApplyCut(wPtCut, "wPtCut"))
        return true;

    // Require OnShell W
    bool onShellW = (recoW.m() >= WMASS_LOW && recoW.m() <= WMASS_HIGH);
    if (!Manager()->ApplyCut(onShellW, "onShellW"))
        return true;

    // Btag Veto
    bool noBTag = BTagVeto(Ak4jets);
    if (not Manager()->ApplyCut(noBTag, "noBTag"))
        return true;

    // dPhi(W, SUEP)
    bool dPhi_W_SUEP = dPhiCut(recoW.phi(), Ak15Jets.at(0).phi(), MIN_DPHI);
    if (not Manager()->ApplyCut(dPhi_W_SUEP, "dPhi_W_SUEP"))
        return true;

    // dPhi(MET, SUEP)
    bool dPhi_MET_SUEP = dPhiCut(event.rec()->MET().phi(), Ak15Jets.at(0).phi(), MIN_DPHI);
    if (not Manager()->ApplyCut(dPhi_MET_SUEP, "dPhi_MET_SUEP"))
        return true;

    // dPhi(Lepton, SUEP)
    bool dPhi_Lep_SUEP = dPhiCut(tight_leptons.at(0).phi(), Ak15Jets.at(0).phi(), MIN_DPHI);
    if (not Manager()->ApplyCut(dPhi_Lep_SUEP, "dPhi_Lep_SUEP"))
        return true;

    // Require Ak4 and Ak15 Overlap
    bool Ak4_Ak15_Overlap = (calculateDeltaR(Ak4jets.at(0), Ak15Jets.at(0)) < 0.4);
    if (not Manager()->ApplyCut(Ak4_Ak15_Overlap, "Ak4_Ak15_Overlap"))
        return true;

    // WpT / SUEPpT < 3
    bool W_SUEP_ratio = (recoW.pt() / Ak15Jets.at(0).pt() < W_SUEP_PT_RATIO);
    if (not Manager()->ApplyCut(W_SUEP_ratio, "W_SUEP_ratio"))
        return true;

    // dPhi(MET, Ak4)
    bool dPhi_MET_Ak4 = dPhiCut(event.rec()->MET().phi(), Ak4jets.at(0).phi(), MIN_DPHI);
    if (not Manager()->ApplyCut(dPhi_MET_Ak4, "dPhi_MET_Ak4"))
        return true;

    //////////////////////////////////////////////
    //           Filling Histograms             //
    //////////////////////////////////////////////

    // Calculating Lab and Boosted Sphericity
    double labSphericity = sphericity(Ak15jetConstituents.at(0), 1.0);
    std::vector<fastjet::PseudoJet> boostedConstituents = boostToSUEP(Ak15result.second.at(0), Ak15result.first.at(0));
    double boostedSphericity = sphericity(boostedConstituents, 1.0);

    // W Histograms
    Manager()->FillHisto("wMass", recoW.m());
    Manager()->FillHisto("wPt", recoW.pt());
    Manager()->FillHisto("wEta", recoW.eta());
    Manager()->FillHisto("wPhi", recoW.phi());

    // Lepton Histograms
    Manager()->FillHisto("lepPt", tight_leptons.at(0).pt());
    Manager()->FillHisto("lepEta", tight_leptons.at(0).eta());
    Manager()->FillHisto("lepPhi", tight_leptons.at(0).phi());
    Manager()->FillHisto("looseLep", loose_leptons.size());
    Manager()->FillHisto("tightLep", tight_leptons.size());

    // Muon or Electron Histograms
    if (tight_muons.size() > 0)
    {
        Manager()->FillHisto("muPt", tight_muons.at(0).pt());
        Manager()->FillHisto("muEta", tight_muons.at(0).eta());
        Manager()->FillHisto("muPhi", tight_muons.at(0).phi());
        Manager()->FillHisto("looseMu", loose_muons.size());
        Manager()->FillHisto("tightMu", tight_muons.size());
    }
    else
    {
        Manager()->FillHisto("elePt", tight_electrons.at(0).pt());
        Manager()->FillHisto("eleEta", tight_electrons.at(0).eta());
        Manager()->FillHisto("elePhi", tight_electrons.at(0).phi());
        Manager()->FillHisto("looseEle", loose_electrons.size());
        Manager()->FillHisto("tightEle", tight_electrons.size());
    }

    // Ak4 Histograms
    Manager()->FillHisto("NJets", Ak4jets.size());
    Manager()->FillHisto("ak4Pt", Ak4jets.at(0).pt());
    Manager()->FillHisto("ak4Eta", Ak4jets.at(0).eta());
    Manager()->FillHisto("ak4Phi", Ak4jets.at(0).phi());
    Manager()->FillHisto("ak4NTracks", Ak4jets.at(0).ntracks());

    // Ak15 Histograms
    Manager()->FillHisto("ak15Pt", Ak15Jets.at(0).pt());
    Manager()->FillHisto("ak15NTracks", Ak15jetConstituents.at(0).size());
    Manager()->FillHisto("ak15Eta", Ak15Jets.at(0).eta());
    double phi_mapped = Ak15Jets.at(0).phi();
    if (phi_mapped > M_PI)
    {
        phi_mapped -= 2 * M_PI;
    }
    Manager()->FillHisto("ak15Phi", phi_mapped);
    Manager()->FillHisto("ak15Mass", Ak15Jets.at(0).m());

    // MET Histograms
    Manager()->FillHisto("metPt", (event.rec)()->MET().pt());
    Manager()->FillHisto("metPhi", event.rec()->MET().phi());
    Manager()->FillHisto("metEta", event.rec()->MET().eta());

    // Sphericity
    Manager()->FillHisto("labSphericity", labSphericity);
    Manager()->FillHisto("boostedSphericity", boostedSphericity);

    // Extended ABCD regions
    if ((10 <= (Ak15jetConstituents.at(0).size()) && (Ak15jetConstituents.at(0).size()) < 20) && (0.3 <= boostedSphericity && boostedSphericity < 0.4))
    {
        Manager()->FillHisto("ABCD_A", Ak15jetConstituents.at(0).size());
    }
    if ((20 <= (Ak15jetConstituents.at(0).size()) && (Ak15jetConstituents.at(0).size()) < 30) && (0.3 <= boostedSphericity && boostedSphericity < 0.4))
    {
        Manager()->FillHisto("ABCD_B", Ak15jetConstituents.at(0).size());
    }
    if ((30 <= (Ak15jetConstituents.at(0).size()) && (Ak15jetConstituents.at(0).size()) < 200) && (0.3 <= boostedSphericity && boostedSphericity < 0.4))
    {
        Manager()->FillHisto("ABCD_C", Ak15jetConstituents.at(0).size());
    }

    if ((10 <= (Ak15jetConstituents.at(0).size()) && (Ak15jetConstituents.at(0).size()) < 20) && (0.4 <= boostedSphericity && boostedSphericity < 0.5))
    {
        Manager()->FillHisto("ABCD_D", Ak15jetConstituents.at(0).size());
    }
    if ((20 <= (Ak15jetConstituents.at(0).size()) && (Ak15jetConstituents.at(0).size()) < 30) && (0.4 <= boostedSphericity && boostedSphericity < 0.5))
    {
        Manager()->FillHisto("ABCD_E", Ak15jetConstituents.at(0).size());
    }
    if ((30 <= (Ak15jetConstituents.at(0).size()) && (Ak15jetConstituents.at(0).size()) < 40) && (0.4 <= boostedSphericity && boostedSphericity < 0.5))
    {
        Manager()->FillHisto("ABCD_F0", Ak15jetConstituents.at(0).size());
    }
    if ((40 <= (Ak15jetConstituents.at(0).size()) && (Ak15jetConstituents.at(0).size()) < 50) && (0.4 <= boostedSphericity && boostedSphericity < 0.5))
    {
        Manager()->FillHisto("ABCD_F1", Ak15jetConstituents.at(0).size());
    }
    if ((50 <= (Ak15jetConstituents.at(0).size()) && (Ak15jetConstituents.at(0).size()) < 60) && (0.4 <= boostedSphericity && boostedSphericity < 0.5))
    {
        Manager()->FillHisto("ABCD_F2", Ak15jetConstituents.at(0).size());
    }
    if ((60 <= (Ak15jetConstituents.at(0).size()) && (Ak15jetConstituents.at(0).size()) < 80) && (0.4 <= boostedSphericity && boostedSphericity < 0.5))
    {
        Manager()->FillHisto("ABCD_F3", Ak15jetConstituents.at(0).size());
    }
    if ((80 <= (Ak15jetConstituents.at(0).size()) && (Ak15jetConstituents.at(0).size()) < 200) && (0.4 <= boostedSphericity && boostedSphericity < 0.5))
    {
        Manager()->FillHisto("ABCD_F4", Ak15jetConstituents.at(0).size());
    }

    if ((10 <= (Ak15jetConstituents.at(0).size()) && (Ak15jetConstituents.at(0).size()) < 20) && (0.5 <= boostedSphericity && boostedSphericity <= 1.0))
    {
        Manager()->FillHisto("ABCD_G", Ak15jetConstituents.at(0).size());
    }
    if ((20 <= (Ak15jetConstituents.at(0).size()) && (Ak15jetConstituents.at(0).size()) < 30) && (0.5 <= boostedSphericity && boostedSphericity <= 1.0))
    {
        Manager()->FillHisto("ABCD_H", Ak15jetConstituents.at(0).size());
    }
    if ((30 <= (Ak15jetConstituents.at(0).size()) && (Ak15jetConstituents.at(0).size()) < 40) && (0.5 <= boostedSphericity && boostedSphericity <= 1.0))
    {
        Manager()->FillHisto("ABCD_SR0", Ak15jetConstituents.at(0).size());
    }
    if ((40 <= (Ak15jetConstituents.at(0).size()) && (Ak15jetConstituents.at(0).size()) < 50) && (0.5 <= boostedSphericity && boostedSphericity <= 1.0))
    {
        Manager()->FillHisto("ABCD_SR1", Ak15jetConstituents.at(0).size());
    }
    if ((50 <= (Ak15jetConstituents.at(0).size()) && (Ak15jetConstituents.at(0).size()) < 60) && (0.5 <= boostedSphericity && boostedSphericity <= 1.0))
    {
        Manager()->FillHisto("ABCD_SR2", Ak15jetConstituents.at(0).size());
    }
    if ((60 <= (Ak15jetConstituents.at(0).size()) && (Ak15jetConstituents.at(0).size()) < 80) && (0.5 <= boostedSphericity && boostedSphericity <= 1.0))
    {
        Manager()->FillHisto("ABCD_SR3", Ak15jetConstituents.at(0).size());
    }
    if ((80 <= (Ak15jetConstituents.at(0).size()) && (Ak15jetConstituents.at(0).size()) < 200) && (0.5 <= boostedSphericity && boostedSphericity <= 1.0))
    {
        Manager()->FillHisto("ABCD_SR4", Ak15jetConstituents.at(0).size());
    }

    return true;
}

///////////////////////////////////////////////////////////////
//                        Finalize                           //
//    Function called one time at the end of the analysis    //
///////////////////////////////////////////////////////////////

void user::Finalize(const SampleFormat &summary, const std::vector<SampleFormat> &files) {}

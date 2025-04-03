#include "SampleAnalyzer/User/Analyzer/user.h"
using namespace MA5;
using namespace std;
#include <Eigen/Dense>
 
// Returns the mass (in GeV) for a charged particle based on its PDG id
// Common charged leptons, mesons, and baryons are included
// If the PDG id is not recognized, defaults to the charged pion mass
double getMassFromPDG(int pdgid) {
    int absid = std::abs(pdgid);
    switch(absid) {
        case 11:    return 0.0005109989461; // electron
        case 13:    return 0.1056583745;    // muon
        case 211:   return 0.13957039;      // charged pion
        case 321:   return 0.493677;        // charged kaon
        case 213:   return 0.77526;         // charged ρ meson
        case 323:   return 0.89166;         // K*(892)+ meson
        case 2212:  return 0.9382720813;    // proton
        case 2214:  return 1.232;           // Δ+ baryon
        case 2224:  return 1.232;           // Δ++ baryon
        case 411:   return 1.86965;         // D+ meson
        case 431:   return 1.96834;         // D_s+ meson
        case 3222:  return 1.18937;         // Σ+ baryon
        case 3112:  return 1.19745;         // Σ- baryon
        case 3312:  return 1.32171;         // Ξ- baryon
        case 3334:  return 1.67245;         // Ω- baryon
        case 521:   return 5.279;           // B+ meson
        case 4122:  return 2.28646;         // Λ_c+ baryon
        case 4222:  return 2.453;           // Σ_c++ baryon
        default: {
            std::cout << "[CAUTION] PDG id " << pdgid 
                      << " of track being clustered into Ak15 not recognized in getMassFromPDG - please update function with this particle's mass to avoid inaccurate results." 
                      << std::endl;
            return 0.13957039;       // default to charged pion mass
        }
    }
}

std::vector<fastjet::PseudoJet> boostToSUEP(const std::vector<fastjet::PseudoJet>& constituents, const fastjet::PseudoJet& jet) 
{
    std::vector<fastjet::PseudoJet> boosted_particles;
    // getting jet four momenta components
    double E_jet = jet.E();
    double px_jet = jet.px();
    double py_jet = jet.py();
    double pz_jet = jet.pz();

    // x, y, and z boosts of the jet
    double bx = px_jet / E_jet;
    double by = py_jet / E_jet;
    double bz = pz_jet / E_jet;

    double beta2 = bx*bx + by*by + bz*bz;
    double gamma = 1.0 / std::sqrt(1.0 - beta2);

    // Debugging Inputs
    //std::cout << "[DEBUG] Jet: E = " << E_jet << ", px = " << px_jet << ", py = " << py_jet << ", pz = " << pz_jet << std::endl;
    //std::cout << "[DEBUG] Boost parameters: bx = " << bx << ", by = " << by << ", bz = " << bz << ", beta2 = " << beta2  << ", gamma = " << gamma << std::endl;

    // boosting constituents
    for (const auto& p : constituents) {
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

        // Debugging Boosting
        //std::cout << "[DEBUG] Dot product (pb) = " << pb << std::endl;
        //std::cout << "[DEBUG] Constituent original: E = " << E << ", px = " << px << ", py = " << py << ", pz = " << pz << std::endl;
        //std::cout << "[DEBUG] Boosted constituent: E' = " << E_prime << ", px' = " << px_prime << ", py' = " << py_prime << ", pz' = " << pz_prime << std::endl;

        fastjet::PseudoJet boosted(px_prime, py_prime, pz_prime, E_prime);
        boosted_particles.push_back(boosted);
    }

    return boosted_particles;
}

double sphericity(const std::vector<fastjet::PseudoJet>& particles, double r) {

    if (particles.empty()) {
        // std::cerr << "Leading jet empty!" << std::endl;
        return 0.0;
    }

    // Initialize sums for sphericity matrix calculation
    double S_xx = 0.0, S_xy = 0.0, S_xz = 0.0;
    double S_yy = 0.0, S_yz = 0.0, S_zz = 0.0;
    double norm = 0.0;

    // Calculate momentum components & and normalization factor
    for (const auto & particle : particles) {
        double px = particle.px();
        double py = particle.py();
        double pz = particle.pz();
        double p = std::sqrt( px * px + py * py + pz * pz);

        if (p == 0) {
            // std::cerr << "Warning: Found particle with zero momentum! Skipping." << std::endl;
            continue;
        }
        
        double weight = std::pow(p, r - 2.0); // Weight for each component

        // Calculating the tensor components
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
    if (norm > 0) {
        S_xx /= norm;
        S_xy /= norm;
        S_xz /= norm;
        S_yy /= norm;
        S_yz /= norm;
        S_zz /= norm;
    }

    // Build sphericity matrix
    Eigen::Matrix3d S;
    S << S_xx, S_xy, S_xz,
         S_xy, S_yy, S_yz,
         S_xz, S_yz, S_zz;

    // Debugging S matrix and normalization
    //std::cout << "[DEBUG] S_xx: " << S_xx << ", S_xy: " << S_xy << ", S_xz: " << S_xz << std::endl;
    //std::cout << "[DEBUG] S_yy: " << S_yy << ", S_yz: " << S_yz << ", S_zz: " << S_zz << std::endl;
    //std::cout << "[DEBUG] Norm: " << norm << std::endl;
    //std::cout << "[DEBUG] Constructed S matrix:\n" << S << std::endl;

    // Get sphericity matrix eigenvalues
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> solver(S);
    if (solver.info() != Eigen::Success) {
        //std::cerr << "[DEBUG] Sphericity eigenvalue solver failed with info: " << solver.info() << std::endl;
        return 0.0;
    }

    // Sort the eigenvalues
    Eigen::Vector3d eigenvals = solver.eigenvalues();
    std::sort(eigenvals.data(), eigenvals.data() + eigenvals.size());

    // Grab the two smallest eigenvalues

    double eigenval1 = eigenvals[0];
    double eigenval2 = eigenvals[1];

    // Calculate the sphericity

    double sphericityval = 1.5 * (eigenval1 + eigenval2);

    // Debugging statements
    // std::cout << "Sphericity matrix: " << S << std::endl;
    // std::cout << "Eigenvalues: " << eigenvals.transpose() << std::endl;
    // std::cout << "Sorted eigenvalues: " << eigenvals[0] << ", " << eigenvals[1] << ", " << eigenvals[2] << std::endl;

    return sphericityval;
}

bool iso(const RecLeptonFormat& lepton, 
    const std::vector<RecTrackFormat>& eflowTracks, 
    const std::vector<RecParticleFormat>& eflowPhotons, 
    const std::vector<RecParticleFormat>& eflowNeutralHadrons, 
    double iso_minpt, 
    double deltaRmax, 
    std::string iso = ""){
   
    // Iso Debugging
    // std::cout << "Inside iso()" << std::endl;

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

    // Iso Debugging
    // std::cout << "ratiomax has been set to: " << ratiomax << std::endl;
    // std::cout << "Tracks size: " << tracks.size() << std::endl;
    // std::cout << "Photons size: " << eflowPhotons.size() << std::endl;
    // std::cout << "Neutral Hadrons size: " << eflowNeutralHadrons.size() << std::endl;

    double lep_pt = lepton.pt();
    double totalpt = 0.0;

    for (const auto & track : eflowTracks) {
        double dr = lepton.dr(track);
        double pt = track.pt();
        //std::cout << "[Track] pt: " << pt << ", dR: " << dr << std::endl;
        if (dr < 0.001) {
            //std::cout << " --> Skipping self" << std::endl;
            continue;
        }
        if (dr < deltaRmax && pt > iso_minpt) {
            //std::cout << " --> Added to totalpt" << std::endl;
            totalpt += pt;
        }
    }
    
    for (const auto & photon : eflowPhotons) {
        double dr = lepton.dr(photon);
        double pt = photon.pt();
        //std::cout << "[Photon] pt: " << pt << ", dR: " << dr << std::endl;
        if (dr < deltaRmax && pt > iso_minpt) {
            //std::cout << " --> Added to totalpt" << std::endl;
            totalpt += pt;
        }
    }
    
    for (const auto & neutral : eflowNeutralHadrons) {
        double dr = lepton.dr(neutral);
        double pt = neutral.pt();
        //std::cout << "[NeutralHadron] pt: " << pt << ", dR: " << dr << std::endl;
        if (dr < deltaRmax && pt > iso_minpt) {
            // std::cout << " --> Added to totalpt" << std::endl;
            totalpt += pt;
        }
    }
    
    double pt_ratio = totalpt / lep_pt;
    //std::cout << "Lepton pt: " << lep_pt << ", total iso pt: " << totalpt << ", pt_ratio: " << pt_ratio << ", ratiomax: " << ratiomax << std::endl;
    return (pt_ratio <= ratiomax);
}

vector<RecLeptonFormat> filter_muons(const vector<RecLeptonFormat>& objects, 
    float ptmin, 
    float etamax, 
    float d0, 
    float dz, 
    const vector<RecTrackFormat>& eflowTracks, 
    const vector<RecParticleFormat>& eflowPhotons, 
    const vector<RecParticleFormat>& eflowNeutralHadrons, 
    const string& selection, 
    string charge = "")
{
    vector<RecLeptonFormat> filtered;
    for (const auto & obj : objects) 
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

        // pT cut
        if (obj.pt() < ptmin)
            continue;

        // For loose selection, apply additional eta and d0 cuts
        if (selection == "loose") 
        {
            if (fabs(obj.eta()) > etamax)
                continue;

            if (fabs(obj.d0()) > d0)
                continue;
        }

        // dz cut applies for both loose and tight selections
        if (fabs(obj.dz()) > dz)
            continue;

        // Iso selection: use "pfIso2" for loose and "pfIso5" for tight
        string iso_option = (selection == "loose") ? "pfIso2" : "pfIso5";
        if (!iso(obj, eflowTracks, eflowPhotons, eflowNeutralHadrons, 0.1, 0.4, iso_option))
            continue;

        filtered.push_back(obj);
    }
    return filtered;
}


vector<RecLeptonFormat> filter_electrons(const vector<RecLeptonFormat>& objects, 
    float ptmin, 
    float etamax, 
    const vector<RecTrackFormat>& eflowTracks, 
    const vector<RecParticleFormat>& eflowPhotons, 
    const vector<RecParticleFormat>& eflowNeutralHadrons, 
    const string& selection, 
    string charge = "")
{
    vector<RecLeptonFormat> filtered;
    for (const auto & obj : objects) 
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

        // pT cut
        if (obj.pt() < ptmin)
        {
            continue;    
        }

        if (selection == "loose") 
        {
            double abseta = fabs(obj.eta());
            if (abseta > etamax)
            {
                continue;
            }
                

            // Exclude crack region: electrons with 1.444 < |eta| < 1.566
            if ((abseta > 1.444) && (abseta < 1.566))
            {
                continue;
            }

            if (fabs(obj.d0()) > (0.05 + 0.05 * (abseta > 1.479)))
            {
                continue;
            }

            if (fabs(obj.dz()) > (0.1 + 0.1 * (abseta > 1.479)))
            {
                continue;
            }

            if (!iso(obj, eflowTracks, eflowPhotons, eflowNeutralHadrons, 0.1, 0.3, "WP90"))
            {
                continue;
            }
        } 

        else if (selection == "tight") 
        {
            if (!iso(obj, eflowTracks, eflowPhotons, eflowNeutralHadrons, 0.1, 0.3, "WP80"))
            {
                continue;
            }
        }

        filtered.push_back(obj);
    }
    return filtered;
}

vector<RecJetFormat> filtered_jets(vector<RecJetFormat> jets, 
    float ptmin, 
    float etamax, 
    float deltaR,
    const vector<RecLeptonFormat>& leptons) {
    // Selects jets that pass analysis-level cuts

    vector<RecJetFormat> filtered;

    for ( auto & jet : jets) {

        // Object Selections

        // pT
        if (jet.pt() < ptmin) {
            continue;
        }

        // eta
        if (fabs(jet.eta()) > etamax) {
            continue;
        }

        // deltaR
        bool passDeltaR = true;

        for (const auto & lepton : leptons) {
            if (jet.dr(lepton) <= deltaR) {
                passDeltaR = false;
                break;
            }
        }

        // putting successful jets in the filtered list
        if (passDeltaR) {
            filtered.push_back(jet);
        }

    }

    return filtered;
        
}

// Orders leptons in leading pT, 0th is leading lepton

bool leadinglep(const RecLeptonFormat &a, const RecLeptonFormat &b) {
    return a.pt() > b.pt();
}

// Orders jets in leading pT, 0th is leading jet

bool leadingjet(const RecJetFormat &a, const RecJetFormat &b) {
    return a.pt() > b.pt();
}

// Sorting the leptons & jets

void leading(std::vector<RecLeptonFormat>& electrons,
             std::vector<RecLeptonFormat>& posElectrons,
             std::vector<RecLeptonFormat>& negElectrons,
             std::vector<RecLeptonFormat>& muons,
             std::vector<RecLeptonFormat>& posMuons,
             std::vector<RecLeptonFormat>& negMuons,
             std::vector<RecLeptonFormat>& leptons,
             std::vector<RecLeptonFormat>& posLeptons,
             std::vector<RecLeptonFormat>& negLeptons,
             std::vector<RecJetFormat> Ak4jets = {}) {
            
    // Sorting electrons
    std::sort(electrons.begin(), electrons.end(), leadinglep);
    std::sort(posElectrons.begin(), posElectrons.end(), leadinglep);
    std::sort(negElectrons.begin(), negElectrons.end(), leadinglep);

    // Sorting muons
    std::sort(muons.begin(), muons.end(), leadinglep);
    std::sort(posMuons.begin(), posMuons.end(), leadinglep);
    std::sort(negMuons.begin(), negMuons.end(), leadinglep);

    // Sorting leptons
    std::sort(leptons.begin(), leptons.end(), leadinglep);
    std::sort(posLeptons.begin(), posLeptons.end(), leadinglep);
    std::sort(negLeptons.begin(), negLeptons.end(), leadinglep);

    // Sorting Ak4 jets
    if (!Ak4jets.empty()) {
        std::sort(Ak4jets.begin(), Ak4jets.end(), leadingjet);
    }
}

// Veto b-tagged jets

bool noBTag(vector<RecJetFormat> jets) {
    for (auto jet : jets) {
        if (jet.btag()) {
            return false;
        }
    }
    return true;
}

// creating the ak15 jet which models the suep

std::pair<std::vector<fastjet::PseudoJet>, std::vector<std::vector<fastjet::PseudoJet>>> getAk15Jets(std::vector <RecTrackFormat> tracks, std::vector<RecLeptonFormat> leptons, double pt_cut, double eta_cut, double dz_cut, double d0_cut, double dr_cut)
{
    std::vector<fastjet::PseudoJet> input_particles;

    for (const auto &track : tracks) {
        if (track.pt() >= pt_cut && std::abs(track.eta()) <= eta_cut &&
            std::abs(track.d0()) < d0_cut && std::abs(track.dz()) < dz_cut && 
            track.dr(leptons.at(0)) >= dr_cut
            )

            // Creating the particles (constituents) which makeup our jet
            {
            double px = track.px();
            double py = track.py();
            double pz = track.pz();
            double p2 = px*px + py*py + pz*pz;
            double mass = getMassFromPDG(track.pdgid());
            double E_corrected = std::sqrt(p2 + mass*mass);
            fastjet::PseudoJet particle(px, py, pz, E_corrected); // Correcting energy of tracks to account for momentum smearing
            input_particles.emplace_back(particle);
            }
    }

    // filling the r = 1.5 jet with its constituents

    fastjet::JetDefinition jet_def(fastjet::antikt_algorithm, 1.5);
    fastjet::ClusterSequence cluster_seq(input_particles, jet_def);

    // again, sort by pt

    std::vector<fastjet::PseudoJet> inclusive_jets = sorted_by_pt(cluster_seq.inclusive_jets(0.0));

    // create a 2d vector which stores jets and their constiuents

    std::vector<std::vector<fastjet::PseudoJet>> cluster_constituents;

    // now we make the jets

    for (const auto &jet : inclusive_jets) 
    {
        cluster_constituents.push_back(jet.constituents());
    }

    return std::make_pair(inclusive_jets, cluster_constituents);
}

// calculating deltaR between the reconstructed jet and the simulated jet

double calculateDeltaR(const RecJetFormat& rec_jet, const fastjet::PseudoJet& pseudo_jet) 
{
    double rec_eta = rec_jet.eta();
    double rec_phi = rec_jet.phi();

    double pseudo_eta = pseudo_jet.eta();
    double pseudo_phi = pseudo_jet.phi();

    // Map pseudo_phi from [0, 2pi] to [-pi, pi] for consistency
    if (pseudo_phi > M_PI) {
        pseudo_phi -= 2 * M_PI;
    }

    // Calculate delta eta and delta phi
    double delta_eta = rec_eta - pseudo_eta;
    double delta_phi = std::atan2(std::sin(rec_phi - pseudo_phi), std::cos(rec_phi - pseudo_phi));

    return std::sqrt(delta_eta * delta_eta + delta_phi * delta_phi);
}

///////////////////////////////////////////////////////////////
//                          init                             //
// function called one time at the beginning of the analysis //
///////////////////////////////////////////////////////////////

MAbool user::Initialize(const MA5::Configuration& cfg, const std::map<std::string, std::string>& parameters)
{
    // initialize the physics service for reconstruction

    PHYSICS->recConfig().Reset();
    PHYSICS->recConfig().UseDeltaRIsolation(0.5);

    // SR definition

    Manager()->AddRegionSelection("SR");

    // Event-Level Selections

    Manager()->AddCut("orthogonality");
    Manager()->AddCut("oneLep");
    Manager()->AddCut("oneAk4");
    Manager()->AddCut("METpTCut_20");
    Manager()->AddCut("Ak15");
    Manager()->AddCut("METpTCut_30");
    Manager()->AddCut("wPtCut");
    Manager()->AddCut("onShellW");
    Manager()->AddCut("noBs");
    Manager()->AddCut("dPhi_W_SUEP");
    Manager()->AddCut("dPhi_MET_SUEP");
    Manager()->AddCut("dPhi_Lep_SUEP");
    Manager()->AddCut("Ak4_Ak15_Overlap");
    Manager()->AddCut("W_SUEP_ratio");
    Manager()->AddCut("dPhi_MET_Ak4");

    // Histos to fill

    // W Histograms
    Manager()->AddHisto("wMass", 60, 0.0, 140.0);
    Manager()->AddHisto("wPt", 300, 0, 300);
    Manager()->AddHisto("wEta", 40,-3.14,3.14);
    Manager()->AddHisto("wPhi", 40,-3.14,3.14);

    // Lepton Histograms
    Manager()->AddHisto("lepPt", 300, 0, 300);
    Manager()->AddHisto("lepEta", 40,-3.14,3.14);
    Manager()->AddHisto("lepPhi", 40,-3.14,3.14);
    Manager()->AddHisto("looseLep", 10, 0, 10);
    Manager()->AddHisto("tightLep", 10, 0, 10);

    // Muon Histograms
    Manager()->AddHisto("muPt", 300, 0, 300);
    Manager()->AddHisto("muEta", 40,-3.14,3.14);
    Manager()->AddHisto("muPhi", 40,-3.14,3.14);
    Manager()->AddHisto("looseMu", 10, 0, 10);
    Manager()->AddHisto("tightMu", 10, 0, 10);

    // Electron Histograms
    Manager()->AddHisto("elePt", 300, 0, 300);
    Manager()->AddHisto("eleEta", 40,-3.14,3.14);
    Manager()->AddHisto("elePhi", 40,-3.14,3.14);
    Manager()->AddHisto("looseEle", 10, 0, 10);
    Manager()->AddHisto("tightEle", 10, 0, 10);

    // Ak4 Histograms
    Manager()->AddHisto("NJets", 10,0.0,10.0);
    Manager()->AddHisto("ak4Pt", 300, 0, 300);
    Manager()->AddHisto("ak4Eta", 40,-3.14,3.14);
    Manager()->AddHisto("ak4Phi", 40,-3.14,3.14);
    Manager()->AddHisto("ak4NTracks", 100,0.0,100.0);

    // Ak15 Histograms
    Manager()->AddHisto("ak15Pt", 300, 0, 300);
    Manager()->AddHisto("ak15Eta", 40,-3.14,3.14);
    Manager()->AddHisto("ak15Phi", 40,-3.14,3.14);
    Manager()->AddHisto("ak15NTracks", 100,0.0,100.0);
    Manager()->AddHisto("ak15Mass", 400,0,400);

    // Sphericity Histograms
    Manager()->AddHisto("labSphericity", 50, 0.0, 1.0);
    Manager()->AddHisto("boostedSphericity", 50, 0.0, 1.0);

    // MET Histograms
    Manager()->AddHisto("metPt", 300, 0, 300);
    Manager()->AddHisto("metPhi", 40, -3.14, 3.14);
    Manager()->AddHisto("metEta", 40, -3.14,3.14);

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
//                          execute                          //
//        function called each time an event is read         //
///////////////////////////////////////////////////////////////

bool user::Execute(SampleFormat& sample, const EventFormat& event)
{   
    // Debugging Event Loop
    //std::cout << "Entering new event" << std::endl;
    //std::cout << "Number of electrons in event: " << event.rec()->electrons().size() << std::endl;
    //std::cout << "Number of muons in event: " << event.rec()->muons().size() << std::endl;
    //std::cout << "Defining object cuts" << std::endl;

    // DEFINING OBJECT CUTS

    // ELECTRONS
    float const LOOSE_ELECTRON_MINPT = 15;
    float const TIGHT_ELECTRON_MINPT = 35;
    float const ELECTRON_MAXETA = 2.5;

    // ELECTRONS NOTES:
    // missing: id, and isolation
    // for delphes card:
        // ####################
        // # Electron isolation
        // ####################

        // module Isolation ElectronIsolation {
        // set CandidateInputArray ElectronEfficiency/electrons
        // set IsolationInputArray EFlowFilter/eflow

        // set OutputArray electrons

        // set DeltaRMax 0.3

        // set PTMin 0.5

        // set PTRatioMax 0.12
        // }

    // MUONS
    float const LOOSE_MUON_MINPT = 10;
    float const TIGHT_MUON_MINPT = 30;
    float const MUON_MAXETA = 2.4;
    float const LOOSE_MUON_DZ = 0.1;
    float const TIGHT_MUON_DZ = 0.05;
    float const MUON_D0 = 0.02;

    // MUONS NOTES:
    // utils ref: https://github.com/SUEPPhysics/SUEPCoffea_dask/blob/d21ebb8d27d8a9c4f264fd680fc2003752aeab8b/workflows/WH_utils.py#L273-L274
    // how to implement: https://github.com/MadAnalysis/madanalysis5/blob/3900ec9002ee7c6963162e72feea58e6971c61eb/tools/SampleAnalyzer/Interfaces/delphes/delphes_cms.tcl#L544

    // Ak4jets
    float const AK4_MINPT = 30;
    float const AK4_MAXETA = 2.4;
    float const AK4LEP_DR = 0.4;

    // Ak4jets NOTES:
    // jet id: https://github.com/SUEPPhysics/SUEPCoffea_dask/blob/d21ebb8d27d8a9c4f264fd680fc2003752aeab8b/workflows/WH_utils.py#L67C73-L67C78

    // APPLYING EVENT WEIGHTS
    // std::cout << "adding weights" << std::endl;
    double weight = 1.;
    if (!Configuration().IsNoEventWeight() && event.mc()!=0) {
        weight = event.mc()->weight();
    }

    Manager()->InitializeForNewEvent(weight);
    if (event.rec() == 0) {return true;}

    //////////////////////////////////////////////
    //  Applying Base Lepton Object Selections  //
    //////////////////////////////////////////////

    // Isolation Eflow Collections
    std::vector<RecTrackFormat> eflowTracks = event.rec()->EFlowTracks();
    std::vector<RecParticleFormat> eflowPhotons = event.rec()->EFlowPhotons();
    std::vector<RecParticleFormat> eflowNeutralHadrons = event.rec()->EFlowNeutralHadrons();

    // Electron Collections

    // Electron Debugging
    //std::cout << "Making electrons" << std::endl;

    vector<RecLeptonFormat> loose_electrons = filter_electrons(event.rec()->electrons(), LOOSE_ELECTRON_MINPT, ELECTRON_MAXETA, eflowTracks, eflowPhotons, eflowNeutralHadrons, "loose");
    vector<RecLeptonFormat> loose_posElectrons = filter_electrons(loose_electrons, LOOSE_ELECTRON_MINPT, ELECTRON_MAXETA, eflowTracks, eflowPhotons, eflowNeutralHadrons, "loose", "+");
    vector<RecLeptonFormat> loose_negElectrons = filter_electrons(loose_electrons, LOOSE_ELECTRON_MINPT, ELECTRON_MAXETA, eflowTracks, eflowPhotons, eflowNeutralHadrons, "loose", "-");

    vector<RecLeptonFormat> tight_electrons = filter_electrons(loose_electrons, TIGHT_ELECTRON_MINPT, ELECTRON_MAXETA, eflowTracks, eflowPhotons, eflowNeutralHadrons, "tight");
    vector<RecLeptonFormat> tight_posElectrons = filter_electrons(tight_electrons, TIGHT_ELECTRON_MINPT, ELECTRON_MAXETA, eflowTracks, eflowPhotons, eflowNeutralHadrons, "tight", "+");
    vector<RecLeptonFormat> tight_negElectrons = filter_electrons(tight_electrons, TIGHT_ELECTRON_MINPT, ELECTRON_MAXETA, eflowTracks, eflowPhotons, eflowNeutralHadrons, "tight", "-");
    
    // Electron Debugging
    //std::cout << "Original electrons size: " << event.rec()->electrons().size() << std::endl;
    //std::cout << "loose_electrons size: " << loose_electrons.size() << std::endl;
    //std::cout << "tight_electrons size: " << tight_electrons.size() << std::endl;

    // Muon Collections

    // Muon Debugging
    //std::cout << "Making muons" << std::endl;

    vector<RecLeptonFormat> loose_muons = filter_muons(event.rec()->muons(), LOOSE_MUON_MINPT, MUON_MAXETA, MUON_D0, LOOSE_MUON_DZ,  eflowTracks, eflowPhotons, eflowNeutralHadrons, "loose");
    vector<RecLeptonFormat> loose_posMuons = filter_muons(loose_muons, LOOSE_MUON_MINPT, MUON_MAXETA, MUON_D0, LOOSE_MUON_DZ, eflowTracks, eflowPhotons, eflowNeutralHadrons,  "loose", "+");
    vector<RecLeptonFormat> loose_negMuons = filter_muons(loose_muons, LOOSE_MUON_MINPT, MUON_MAXETA, MUON_D0, LOOSE_MUON_DZ, eflowTracks, eflowPhotons, eflowNeutralHadrons,  "loose", "-");

    vector<RecLeptonFormat> tight_muons = filter_muons(loose_muons, TIGHT_MUON_MINPT, MUON_MAXETA, MUON_D0, TIGHT_MUON_DZ, eflowTracks, eflowPhotons, eflowNeutralHadrons, "tight");
    vector<RecLeptonFormat> tight_posMuons = filter_muons(tight_muons, TIGHT_MUON_MINPT, MUON_MAXETA, MUON_D0, TIGHT_MUON_DZ, eflowTracks, eflowPhotons, eflowNeutralHadrons, "tight", "+");
    vector<RecLeptonFormat> tight_negMuons = filter_muons(tight_muons, TIGHT_MUON_MINPT, MUON_MAXETA, MUON_D0, TIGHT_MUON_DZ, eflowTracks, eflowPhotons, eflowNeutralHadrons, "tight", "-");

    // Muon Debugging
    //std::cout << "Original muons size: " << event.rec()->muons().size() << std::endl;
    //std::cout << "loose_muons size: " << loose_muons.size() << std::endl;
    //std::cout << "tight_muons size: " << tight_muons.size() << std::endl;

    // Combining into loose Leptons
    vector<RecLeptonFormat> loose_leptons = loose_electrons;
    loose_leptons.insert(loose_leptons.end(), loose_muons.begin(), loose_muons.end());

    vector<RecLeptonFormat> loose_posLeptons = loose_posElectrons;
    loose_posLeptons.insert(loose_posLeptons.end(), loose_posMuons.begin(), loose_posMuons.end());

    vector<RecLeptonFormat> loose_negLeptons = loose_negElectrons;
    loose_negLeptons.insert(loose_negLeptons.end(), loose_negMuons.begin(), loose_negMuons.end());

    // Combining into tight Leptons
    vector<RecLeptonFormat> tight_leptons = tight_electrons;
    tight_leptons.insert(tight_leptons.end(), tight_muons.begin(), tight_muons.end());

    vector<RecLeptonFormat> tight_posLeptons = tight_posElectrons;
    tight_posLeptons.insert(tight_posLeptons.end(), tight_posMuons.begin(), tight_posMuons.end());

    vector<RecLeptonFormat> tight_negLeptons = tight_negElectrons;
    tight_negLeptons.insert(tight_negLeptons.end(), tight_negMuons.begin(), tight_negMuons.end());

    // Lepton Debugging
    //std::cout << "Loose Leptons size: " << loose_leptons.size() << std::endl;
    //std::cout << "Tight Leptons size: " << tight_leptons.size() << std::endl;

    // Sorting leptons
    leading(loose_electrons, loose_posElectrons, loose_negElectrons, loose_muons, loose_posMuons, loose_negMuons, loose_leptons, loose_posLeptons, loose_negLeptons);
    leading(tight_electrons, tight_posElectrons, tight_negElectrons, tight_muons, tight_posMuons, tight_negMuons, tight_leptons, tight_posLeptons, tight_negLeptons);

    //////////////////////////////////////////////
    //           Event level selections         //
    //////////////////////////////////////////////

    //std::cout << "Starting event level selections" << std::endl;

    // Orthogonality to GGF offline: Remove events with no leptons
    //std::cout << "GGF Orthogonality" << std::endl;
    bool GGFOrthogonality = (loose_leptons.size() == 0);

    // Orthogonality to ZH: Remove events with a pair of loose OSSF leptons
    // std::cout << "ZH Orthogonality" << std::endl;
    bool twoOSleptons = (loose_posLeptons.size() == 1 && loose_negLeptons.size() == 1);                  // Require exactly two opposite-sign leptons (either muons or electrons)
    bool twoSFleptons = (loose_muons.size() == 2 || loose_electrons.size() == 2);                        // Require exactly two muons or exactly two electrons
    bool LeadpTleptons = (loose_leptons.size() > 0) && (loose_leptons.at(0).pt() >= 25);                 // Require the leading lepton pT to be >= 25 GeV
    bool ZHOrthogonality = twoOSleptons && twoSFleptons && LeadpTleptons;                                // Concatenating cuts
    
    // Apply both orthogonality cuts
    //std::cout << "Joint orthogonality" << std::endl;
    bool orthogonality = (!GGFOrthogonality) && (!ZHOrthogonality);
    if (not Manager()->ApplyCut(orthogonality, "orthogonality")) return true;

    // Exactly one tight lepton
    //std::cout << "Exactly One Tight" << std::endl;
    bool oneLep = (tight_leptons.size() == 1);
    if (not Manager()->ApplyCut(oneLep, "oneLep")) return true;

    // Debugging Leptons
    //std::cout << "Tight Electrons size: " << tight_electrons.size() << std::endl;
    //std::cout << "Tight Muons size: " << tight_muons.size() << std::endl;
    //std::cout << "Tight Leptons size: " << tight_leptons.size() << std::endl;
    // std::cout << "Lepton px: " << tight_leptons.at(0).px() 
    //         << " py: " << tight_leptons.at(0).py() 
    //         << " pz: " << tight_leptons.at(0).pz() 
    //         << " E: " << tight_leptons.at(0).e() << std::endl;
    
    // Ak4 Jet Collection
    //std::cout << "Defining ak4jet objects" << std::endl;
    vector<RecJetFormat> Ak4jets = filtered_jets(event.rec()->jets(), AK4_MINPT, AK4_MAXETA, AK4LEP_DR, tight_leptons);

    // Sorting Leptons and Ak4s
    leading(tight_electrons, tight_posElectrons, tight_negElectrons, tight_muons, tight_posMuons, tight_negMuons, tight_leptons, tight_posLeptons, tight_negLeptons, Ak4jets);

    // Require at least one ak4
    bool oneAk4 = (Ak4jets.size() > 0);
    if (not Manager()->ApplyCut(oneAk4, "oneAk4")) return true;

    // Require MET > 20 GeV
    bool METpTCut_20 = (event.rec()->MET().pt() > 20);
    if (!Manager()->ApplyCut(METpTCut_20, "METpTCut_20")) return true;

    //////////////////////////////////////////////
    //  Building the W and Clustering the Ak15  //
    //////////////////////////////////////////////

    // W Boson
    float const MIN_W_PT = 60;
    float const MIN_W_MASS = 30;
    float const MAX_W_MASS = 130;

    // Track ID
    float const MIN_TRACK_PT = 1;
    float const MAX_TRACK_ETA = 2.5;
    float const MAX_TRACK_D0 = 0.05;
    float const MAX_TRACK_DZ = 0.05;
    float const MIN_TRACK_LEP_DR = 0.4;

    // Ak15 pT
    float const AK15JET_PT_MIN = 60;

    // Misc
    float const MIN_DPHI = 1.5;
    float const W_SUEP_PT_RATIO = 3.0;

    // W reconstruction
    ParticleBaseFormat recoW;
    recoW += tight_leptons.at(0).momentum();
    recoW += event.rec()->MET().momentum();

    // W Debugging
    // std::cout << "W candidate mass: " << recoW.m() << std::endl;

    // Jet Clustering
    //std::cout << "Building Ak15 Jets" << std::endl;
    std::pair<std::vector<fastjet::PseudoJet>, std::vector<std::vector<fastjet::PseudoJet>>> Ak15Result = getAk15Jets(event.rec()->tracks(), tight_leptons, MIN_TRACK_PT, MAX_TRACK_ETA, MAX_TRACK_DZ, MAX_TRACK_D0, MIN_TRACK_LEP_DR);
    //std::cout << "Ak15s clustered" << std::endl;

    std::vector<fastjet::PseudoJet> Ak15Jets = Ak15Result.first;
    std::vector<std::vector<fastjet::PseudoJet>> Ak15jetConstituents = Ak15Result.second;

    // Jet Clustering Debugging
    // std::cout << "Ak15Jets size: " << Ak15Jets.size() << std::endl;
    // std::cout << "Ak15jetConstituents size: " << Ak15jetConstituents.size() << std::endl;

    // std::cout << "Printing Ak15Jets elements:" << std::endl;
    // for (size_t i = 0; i < Ak15Jets.size(); i++) {
    //     std::cout << "Jet " << i << " pt: " << Ak15Jets[i].pt() << std::endl;
    // }

    // std::cout << "Printing Ak15jetConstituents elements:" << std::endl;
    // for (size_t i = 0; i < Ak15jetConstituents.size(); i++) {
    //     if (Ak15jetConstituents[i].empty()) {
    //         std::cerr << "WARNING: Ak15jetConstituents[" << i << "] is empty!" << std::endl;
    //     } else {
    //         std::cout << "Ak15jetConstituents[" << i << "] size: " << Ak15jetConstituents[i].size() << std::endl;
    //     }
    // }


    //////////////////////////////////////////////
    //           Event level selections         //
    //////////////////////////////////////////////

    // Require at least one Ak15 with pT > 60 GeV
    bool oneAk15 = (Ak15Jets.size() > 0);
    bool minAk15pT = (Ak15Jets.at(0).pt() > AK15JET_PT_MIN);
    bool Ak15 = (oneAk15) && (minAk15pT);
    if (not Manager()->ApplyCut(Ak15, "Ak15")) return true;

    // Require MET > 30 GeV
    bool METpTCut_30 = (event.rec()->MET().pt() > 30);
    if (!Manager()->ApplyCut(METpTCut_30, "METpTCut_30")) return true;

    // MET Debugging
    // if (!event.rec()) {
    //     std::cerr << "ERROR: event.rec() is NULL!" << std::endl;
    //     return false;
    // }
    // std::cout << "MET px: " << event.rec()->MET().px() 
    //         << " py: " << event.rec()->MET().py() 
    //         << " pz: " << event.rec()->MET().pz() 
    //         << " E: " << event.rec()->MET().e() << std::endl;

    // Require W pT > 60 GeV
    bool wPtCut = (recoW.pt() > MIN_W_PT);
    if (not Manager()->ApplyCut(wPtCut, "wPtCut")) return true;

    // Require OnShell W
    bool onShellW = (recoW.m() >= MIN_W_MASS && recoW.m() <= MAX_W_MASS);
    if (!Manager()->ApplyCut(onShellW, "onShellW")) return true;

    // Btag Veto
    bool noBs = noBTag(Ak4jets);
    if (not Manager()->ApplyCut(noBs, "noBs")) return true;

    // dPhi(W, SUEP)
    float dPhi_W_SUEP_Val = recoW.phi() - Ak15Jets.at(0).phi();
    if (dPhi_W_SUEP_Val > M_PI) 
    {
        dPhi_W_SUEP_Val = dPhi_W_SUEP_Val - 2*M_PI;
    }
    else if (dPhi_W_SUEP_Val < -M_PI)
    {
        dPhi_W_SUEP_Val = dPhi_W_SUEP_Val + 2*M_PI;
    }
    bool dPhi_W_SUEP = (fabs(dPhi_W_SUEP_Val) > MIN_DPHI);
    if (not Manager()->ApplyCut(dPhi_W_SUEP, "dPhi_W_SUEP")) return true;

    // dPhi(MET, SUEP)
    float dPhi_MET_SUEP_Val = event.rec()->MET().phi() - Ak15Jets.at(0).phi();
    if (dPhi_MET_SUEP_Val > M_PI) 
    {
        dPhi_MET_SUEP_Val = dPhi_MET_SUEP_Val - 2*M_PI;
    }
    else if (dPhi_MET_SUEP_Val < -M_PI)
    {
        dPhi_MET_SUEP_Val = dPhi_MET_SUEP_Val + 2*M_PI;
    }
    bool dPhi_MET_SUEP = (fabs(dPhi_MET_SUEP_Val) > MIN_DPHI);
    if (not Manager()->ApplyCut(dPhi_MET_SUEP, "dPhi_MET_SUEP")) return true;
    
    // dPhi(Lepton, SUEP)
    float dPhi_Lep_SUEP_Val = tight_leptons.at(0).phi() - Ak15Jets.at(0).phi();
    if (dPhi_Lep_SUEP_Val > M_PI) 
    {
        dPhi_Lep_SUEP_Val = dPhi_Lep_SUEP_Val - 2*M_PI;
    }
    else if (dPhi_Lep_SUEP_Val < -M_PI)
    {
        dPhi_Lep_SUEP_Val = dPhi_Lep_SUEP_Val + 2*M_PI;
    }
    bool dPhi_Lep_SUEP = (fabs(dPhi_Lep_SUEP_Val) > MIN_DPHI);
    if (not Manager()->ApplyCut(dPhi_Lep_SUEP, "dPhi_Lep_SUEP")) return true;

    // Require Ak4 and Ak15 Overlap
    bool Ak4_Ak15_Overlap = (calculateDeltaR(Ak4jets.at(0), Ak15Jets.at(0)) < 0.4);
    if (not Manager()->ApplyCut(Ak4_Ak15_Overlap, "Ak4_Ak15_Overlap")) return true;

    // WpT / SUEPpT < 3
    bool W_SUEP_ratio = (recoW.pt() / Ak15Jets.at(0).pt() < W_SUEP_PT_RATIO);
    if (not Manager()->ApplyCut(W_SUEP_ratio, "W_SUEP_ratio")) return true;

    // dPhi(MET, Ak4)
    float dPhi_MET_Ak4_Val = event.rec()->MET().phi() - Ak4jets.at(0).phi();
    if (dPhi_MET_Ak4_Val > M_PI) 
    {
        dPhi_MET_Ak4_Val = dPhi_MET_Ak4_Val - 2*M_PI;
    }
    else if (dPhi_MET_Ak4_Val < -M_PI)
    {
        dPhi_MET_Ak4_Val = dPhi_MET_Ak4_Val + 2*M_PI;
    }
    bool dPhi_MET_Ak4 = (fabs(dPhi_MET_Ak4_Val) > MIN_DPHI);
    if (not Manager()->ApplyCut(dPhi_MET_Ak4, "dPhi_MET_Ak4")) return true;

    // Debugging event selections
    //std::cout << "All cuts successfully applied" << std::endl;

    //////////////////////////////////////////////
    //           Filling Histograms             //
    //////////////////////////////////////////////

    // Calculating Lab and Boosted Sphericity
    double labSphericity = sphericity(Ak15jetConstituents.at(0), 1.0);
    std::vector<fastjet::PseudoJet> boostedConstituents = boostToSUEP(Ak15Result.second.at(0), Ak15Result.first.at(0));
    double boostedSphericity = sphericity(boostedConstituents, 1.0);

    // W Histograms
    Manager()->FillHisto("wMass", recoW.m());
    Manager()->FillHisto("wPt", recoW.pt());
    Manager()->FillHisto("wEta", recoW.eta());
    Manager()->FillHisto("wPhi", recoW.phi());

    // Muons or Electron Histograms
    if(tight_muons.size() > 0) 
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

    // Lepton Histograms
    Manager()->FillHisto("lepPt", tight_leptons.at(0).pt());
    Manager()->FillHisto("lepEta", tight_leptons.at(0).eta());
    Manager()->FillHisto("lepPhi", tight_leptons.at(0).phi());
    Manager()->FillHisto("looseLep", loose_leptons.size());
    Manager()->FillHisto("tightLep", tight_leptons.size());

    // Ak4 Histograms
    Manager()->FillHisto("NJets", Ak4jets.size());
    Manager()->FillHisto("ak4Pt", Ak4jets.at(0).pt());
    Manager()->FillHisto("ak4Eta", Ak4jets.at(0).eta());
    Manager()->FillHisto("ak4Phi", Ak4jets.at(0).phi());
    Manager()->FillHisto("ak4NTracks", Ak4jets.at(0).ntracks());

    // Ak15 Histograms
    Manager()->FillHisto("ak15Pt", Ak15Jets.at(0).pt());
    Manager()->FillHisto("ak15Eta", Ak15Jets.at(0).eta());
    double phi_recalc = Ak15Jets.at(0).phi();
    if (phi_recalc > M_PI){
        phi_recalc -= 2 * M_PI;
    }
    Manager()->FillHisto("ak15Phi", phi_recalc);
    Manager()->FillHisto("ak15NTracks", Ak15jetConstituents.at(0).size());
    Manager()->FillHisto("ak15Mass", Ak15Jets.at(0).m());

    
    // MET Histograms
    Manager()->FillHisto("metPt", (event.rec)()->MET().pt());
    Manager()->FillHisto("metPhi", event.rec()->MET().phi());
    Manager()->FillHisto("metEta", event.rec()->MET().eta());

    // Sphericity
    Manager()->FillHisto("labSphericity", labSphericity);
    Manager()->FillHisto("boostedSphericity", boostedSphericity);

    // extended abcd regions

    if ((10 <= (Ak15jetConstituents.at(0).size()) && (Ak15jetConstituents.at(0).size()) < 20) && (0.3 <= boostedSphericity && boostedSphericity < 0.4)) {Manager()->FillHisto("ABCD_A", Ak15jetConstituents.at(0).size());}
    if ((20 <= (Ak15jetConstituents.at(0).size()) && (Ak15jetConstituents.at(0).size()) < 30) &&  (0.3 <= boostedSphericity && boostedSphericity < 0.4)){Manager()->FillHisto("ABCD_B", Ak15jetConstituents.at(0).size());}
    if ((30 <= (Ak15jetConstituents.at(0).size()) && (Ak15jetConstituents.at(0).size()) < 200) &&  (0.3 <= boostedSphericity && boostedSphericity < 0.4)){Manager()->FillHisto("ABCD_C", Ak15jetConstituents.at(0).size());}

    if ((10 <= (Ak15jetConstituents.at(0).size()) && (Ak15jetConstituents.at(0).size()) < 20) &&  (0.4 <= boostedSphericity && boostedSphericity < 0.5)){Manager()->FillHisto("ABCD_D", Ak15jetConstituents.at(0).size());}
    if ((20 <= (Ak15jetConstituents.at(0).size()) && (Ak15jetConstituents.at(0).size()) < 30) &&  (0.4 <= boostedSphericity && boostedSphericity < 0.5)){Manager()->FillHisto("ABCD_E", Ak15jetConstituents.at(0).size());}
    if ((30 <= (Ak15jetConstituents.at(0).size()) && (Ak15jetConstituents.at(0).size()) < 40) &&  (0.4 <= boostedSphericity && boostedSphericity < 0.5)){Manager()->FillHisto("ABCD_F0", Ak15jetConstituents.at(0).size());}
    if ((40 <= (Ak15jetConstituents.at(0).size()) && (Ak15jetConstituents.at(0).size()) < 50) &&  (0.4 <= boostedSphericity && boostedSphericity < 0.5)){Manager()->FillHisto("ABCD_F1", Ak15jetConstituents.at(0).size());}
    if ((50 <= (Ak15jetConstituents.at(0).size()) && (Ak15jetConstituents.at(0).size()) < 60) &&  (0.4 <= boostedSphericity && boostedSphericity < 0.5)){Manager()->FillHisto("ABCD_F2", Ak15jetConstituents.at(0).size());}
    if ((60 <= (Ak15jetConstituents.at(0).size()) && (Ak15jetConstituents.at(0).size()) < 80) &&  (0.4 <= boostedSphericity && boostedSphericity < 0.5)){Manager()->FillHisto("ABCD_F3", Ak15jetConstituents.at(0).size());}
    if ((80 <= (Ak15jetConstituents.at(0).size()) && (Ak15jetConstituents.at(0).size()) < 200) &&  (0.4 <= boostedSphericity && boostedSphericity < 0.5)){Manager()->FillHisto("ABCD_F4", Ak15jetConstituents.at(0).size());}
    
    if ((10 <= (Ak15jetConstituents.at(0).size()) && (Ak15jetConstituents.at(0).size()) < 20) &&  (0.5 <= boostedSphericity && boostedSphericity <= 1.0)){Manager()->FillHisto("ABCD_G", Ak15jetConstituents.at(0).size());}
    if ((20 <= (Ak15jetConstituents.at(0).size()) && (Ak15jetConstituents.at(0).size()) < 30) &&  (0.5 <= boostedSphericity && boostedSphericity <= 1.0)){Manager()->FillHisto("ABCD_H", Ak15jetConstituents.at(0).size());}
    if ((30 <= (Ak15jetConstituents.at(0).size()) && (Ak15jetConstituents.at(0).size()) < 40) &&  (0.5 <= boostedSphericity && boostedSphericity <= 1.0)){Manager()->FillHisto("ABCD_SR0", Ak15jetConstituents.at(0).size());}
    if ((40 <= (Ak15jetConstituents.at(0).size()) && (Ak15jetConstituents.at(0).size()) < 50) &&  (0.5 <= boostedSphericity && boostedSphericity <= 1.0)){Manager()->FillHisto("ABCD_SR1", Ak15jetConstituents.at(0).size());}
    if ((50 <= (Ak15jetConstituents.at(0).size()) && (Ak15jetConstituents.at(0).size()) < 60) &&  (0.5 <= boostedSphericity && boostedSphericity <= 1.0)){Manager()->FillHisto("ABCD_SR2", Ak15jetConstituents.at(0).size());}
    if ((60 <= (Ak15jetConstituents.at(0).size()) && (Ak15jetConstituents.at(0).size()) < 80) &&  (0.5 <= boostedSphericity && boostedSphericity <= 1.0)){Manager()->FillHisto("ABCD_SR3", Ak15jetConstituents.at(0).size());}
    if ((80 <= (Ak15jetConstituents.at(0).size()) && (Ak15jetConstituents.at(0).size()) < 200) &&  (0.5 <= boostedSphericity && boostedSphericity <= 1.0)){Manager()->FillHisto("ABCD_SR4", Ak15jetConstituents.at(0).size());}

    // std::cout << "All histos successfully filled" << std::endl;
    return true;
}

///////////////////////////////////////////////////////////////
//                        Finalize                           //
//    Function called one time at the end of the analysis    //
///////////////////////////////////////////////////////////////

void user::Finalize(const SampleFormat& summary, const std::vector<SampleFormat>& files){}

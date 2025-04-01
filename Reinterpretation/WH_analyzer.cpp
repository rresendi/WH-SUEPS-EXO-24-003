#include "SampleAnalyzer/User/Analyzer/user.h"
using namespace MA5;
using namespace std;
#include <Eigen/Dense>

std::vector<fastjet::PseudoJet> boost_to_jet_frame(const std::vector<fastjet::PseudoJet>& constituents, const fastjet::PseudoJet& jet) 
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

        fastjet::PseudoJet boosted(px_prime, py_prime, pz_prime, E_prime);
        boosted_particles.push_back(boosted);
    }

    return boosted_particles;
}


double iso(vector<RecLeptonFormat> leptons, std::vector <RecTrackFormat> tracks, std::vector <RecParticleFormat> eflowPhotons, std::vector <RecParticleFormat> eflowNeutralHadrons, double iso_minpt, std::string iso = "")
{
    // std::cout << "Inside iso()" << std::endl;
    // std::cout << "lep size in iso(): " << leptons.size() << std::endl;

    if (leptons.size() == 0) {
        // std::cout << "No leptons to process for iso." << std::endl;
        return 0;
    }

    float deltaRmax = 0.4;
    float ratiomax = 0.0;

    // to be played with
    if (iso == "WP90") {
        ratiomax = 1000000;
        }
    else if (iso == "WP80") {
        ratiomax = 1000000;
        }
    else if (iso == "pfIso2") {
        ratiomax = 1000000;
        }
    else if (iso == "pfIso5") {
        ratiomax = 1000000;
        }

    // std::cout << "ratiomax has been set to: " << ratiomax << std::endl;
    // std::cout << "Tracks size: " << tracks.size() << std::endl;
    // std::cout << "Photons size: " << eflowPhotons.size() << std::endl;
    // std::cout << "Neutral Hadrons size: " << eflowNeutralHadrons.size() << std::endl;


    for (const auto & lepton : leptons) {
        float lep1pt = lepton.pt();
        float totalpt = 0.0;
        for (const auto & track : tracks){
            // std::cout << "track pt: " << track.pt() << "track lep dr: " << lepton.dr(track) << std::endl;
            if ((lepton.dr(track) < deltaRmax) && track.pt() > iso_minpt) {
                // std::cout << "track pt: " << track.pt() << std::endl;
                totalpt += track.pt();
            }
        }
        for (const auto & photon : eflowPhotons){
            // std::cout << "photon pt: " << photon.pt() << "photon lep dr: " << lepton.dr(photon) << std::endl;
            if ((lepton.dr(photon) < deltaRmax) && photon.pt() > iso_minpt) {
                // std::cout << "photon pt: " << photon.pt() << std::endl;
                totalpt += photon.pt();
            }
        }
        for (const auto & neutral : eflowNeutralHadrons){
            // std::cout << "neutral pt: " << neutral.pt() << "neutral lep dr: " << lepton.dr(neutral) << std::endl;
            if ((lepton.dr(neutral) < deltaRmax) && neutral.pt() > iso_minpt) {
                // std::cout << "neutral pt: " << neutral.pt() << std::endl;
                totalpt += neutral.pt();
            }
        }
        // std::cout << "lep1pt: " << lep1pt << ", totalpt: " << totalpt << std::endl;
        float pt_ratio = totalpt / lep1pt;
        // std::cout << "pt_ratio: " << pt_ratio << std::endl;
        if (pt_ratio > ratiomax) {
            // std::cout << "returning 0" << std::endl;
            return 0;
        }
    }
    // std::cout << "returning 1" << std::endl;
    return 1;
}

double sphericity(const std::vector<fastjet::PseudoJet>& particles, double r) {
    // std::cout << "Inside sphericity(). Input vector size: " << particles.size() << ", r = " << r << std::endl;
    // check for an empty jet

    if (particles.empty()) {
        // std::cerr << "Leading jet empty!" << std::endl;
        return 0.0;
    }

    // initialize sums for sphericity matrix calculation

    double S_xx = 0.0, S_xy = 0.0, S_xz = 0.0;
    double S_yy = 0.0, S_yz = 0.0, S_zz = 0.0;
    double norm = 0.0;

    // calculate momentum components & momentum

    for (const auto & particle : particles) {
        double px = particle.px();
        double py = particle.py();
        double pz = particle.pz();
        double p = std::sqrt( px * px + py * py + pz * pz);
        // std::cout << "Particle px: " << px << ", py: " << py << ", pz: " << pz << ", p: " << p << std::endl;

        if (p == 0) {
            // std::cerr << "Warning: Found particle with zero momentum! Skipping." << std::endl;
            continue;
        }
        
        // calculating the weight

        double weight = std::pow(p, r - 2.0);

        // calculating the tensor components

        S_xx += px * px * weight;
        S_xy += px * py * weight;
        S_xz += px * pz * weight;
        S_yy += py * py * weight;
        S_yz += py * pz * weight;
        S_zz += pz * pz * weight;

        // calculating the normalization for the components
        norm += std::pow(p, r);
        // std::cout << "Final norm before normalization: " << norm << std::endl;

    }

    // normalize the matrix if norm > 0
    
    if (norm > 0) {
        S_xx /= norm;
        S_xy /= norm;
        S_xz /= norm;
        S_yy /= norm;
        S_yz /= norm;
        S_zz /= norm;
    }

    // build sphericity matrix

    Eigen::Matrix3d S;
    S << S_xx, S_xy, S_xz,
         S_xy, S_yy, S_yz,
         S_xz, S_yz, S_zz;

    // get sphericity matrix eigenvalues
    // std::cout << "Sphericity matrix: " << S << std::endl;
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> solver(S);
    if (solver.info()!= Eigen::Success) {
        // std::cerr << "Could not calculate eigenvalues..." << std::endl;
        return 0.0;
    }

    // sort the eigenvalues
    Eigen::Vector3d eigenvals = solver.eigenvalues();
    // std::cout << "Eigenvalues: " << eigenvals.transpose() << std::endl;
    std::sort(eigenvals.data(), eigenvals.data() + eigenvals.size());
    // std::cout << "Sorted eigenvalues: " << eigenvals[0] << ", " << eigenvals[1] << ", " << eigenvals[2] << std::endl;

    // grab the two smallest eigenvalues

    double eigenval1 = eigenvals[0];
    double eigenval2 = eigenvals[1];

    // calculate the sphericity

    double sphericityval = 1.5 * (eigenval1 + eigenval2);

    return sphericityval;
}

vector<RecLeptonFormat> filtered_loose_muons(vector<RecLeptonFormat> objects, float ptmin, float etamax, float dz, float d0, std::vector <RecTrackFormat> tracks, std::vector <RecParticleFormat> eflowPhotons, std::vector <RecParticleFormat> eflowNeutralHadrons, std::string charge = "") {
    // Selects loose muons that pass analysis-level cuts

    vector<RecLeptonFormat> filtered;
    // std::cout << "object size: " << objects.size() << std::endl;


    for (auto & obj : objects) {
        
        // Charge Selections
        if ( charge == '+') {
            if (obj.charge() <= 0) {
                continue;
            }
        } else if ( charge == '-') {
            if (obj.charge() >= 0) {
                continue;
            }
        }

        // Object Selections

        // pT
        if (obj.pt() < ptmin) {
            continue;
        }


        // eta
        if (fabs(obj.eta()) > etamax) {
            continue;
        }


        // dz

        if (fabs(obj.dz()) > dz) {
            continue;
        }


        // d0

        if (fabs(obj.d0()) > d0) {
            continue;
        }

        // std::cout << "Before iso call, obj pt: " << obj.pt() << ", eta: " << obj.eta() << std::endl;
        // std::cout << "Tracks size: " << tracks.size() << std::endl;
        // std::cout << "Photons size: " << eflowPhotons.size() << std::endl;
        // std::cout << "Neutral Hadrons size: " << eflowNeutralHadrons.size() << std::endl;

        double iso_value = iso({obj}, tracks, eflowPhotons, eflowNeutralHadrons, 0.1, "pfIso2");
        // std::cout << "iso returned: " << iso_value << std::endl;


        if (iso_value == 0) {
            continue;
        }

        // std::cout << "made it past the iso function in filtered_loose_muons" << std::endl;


        filtered.push_back(obj);
        // std::cout << "successfully appended" << std::endl;


    }
    // std::cout << "returning filtered" << std::endl;
    return filtered;
        
}

vector<RecLeptonFormat> filtered_tight_muons(vector<RecLeptonFormat> objects, float ptmin, float dz, std::vector <RecTrackFormat> tracks, std::vector <RecParticleFormat> eflowPhotons, std::vector <RecParticleFormat> eflowNeutralHadrons, std::string charge = "") {
    // Selects tight muons that pass analysis-level cuts

    vector<RecLeptonFormat> filtered;

    for ( auto & obj : objects) {
        
        // Charge Selections
        if ( charge == '+') {
            if (obj.charge() <= 0) {
                continue;
            }
        } else if ( charge == '-') {
            if (obj.charge() >= 0) {
                continue;
            }
        }

        // Object Selections

        // pT
        if (obj.pt() < ptmin) {
            continue;
        }

        // dz

        if (fabs(obj.dz()) > dz) {
            continue;
        }

        if (iso({obj}, tracks, eflowPhotons, eflowNeutralHadrons, 0.1, "pfIso5") == 0) {
            continue;
        }

        filtered.push_back(obj);

    }

    return filtered;
        
}

vector<RecLeptonFormat> filtered_loose_electrons(vector<RecLeptonFormat> objects, float ptmin, float etamax, std::vector <RecTrackFormat> tracks, std::vector <RecParticleFormat> eflowPhotons, std::vector <RecParticleFormat> eflowNeutralHadrons, std::string charge = "") {
    // Selects loose electrons that pass analysis-level cuts

    vector<RecLeptonFormat> filtered;

    for ( auto & obj : objects) {
        
        // Charge Selections
        if ( charge == '+') {
            if (obj.charge() <= 0) {
                continue;
            }
        } else if ( charge == '-') {
            if (obj.charge() >= 0) {
                continue;
            }
        }

        // Object Selections
        double abseta = fabs(obj.eta());

        // pT
        if (obj.pt() < ptmin) {
            continue;
        }

        // eta and impact parameter
        if (abseta > etamax) {
            continue;
        }

        if (abseta > 1.444 || abseta < 1.566) {
          continue;
        }

        if(fabs(obj.d0()) > (0.05 + 0.05 * (abseta > 1.479))){
          continue;
        }
    
        if(fabs(obj.dz()) > (0.1 + 0.1 * (abseta > 1.479))){
          continue;
        }

        // wp 90

        if (iso({obj}, tracks, eflowPhotons, eflowNeutralHadrons, 0.1, "WP90") == 0) {
            continue;
        }

        filtered.push_back(obj);

    }

    return filtered;
        
}

vector<RecLeptonFormat> filtered_tight_electrons(vector<RecLeptonFormat> objects, float ptmin, std::vector <RecTrackFormat> tracks, std::vector <RecParticleFormat> eflowPhotons, std::vector <RecParticleFormat> eflowNeutralHadrons, std::string charge = "") {
    // Selects tight electrons that pass analysis-level cuts

    vector<RecLeptonFormat> filtered;

    for ( auto & obj : objects) {
        
        // Charge Selections
        if ( charge == '+') {
            if (obj.charge() <= 0) {
                continue;
            }
        } else if ( charge == '-') {
            if (obj.charge() >= 0) {
                continue;
            }
        }

        // Object Selections
        
        // pT
        if (obj.pt() < ptmin) {
            continue;
        }

        // wp 80
        
        if (iso({obj}, tracks, eflowPhotons, eflowNeutralHadrons, 0.1, "WP80") == 0) {
            continue;
        }

        filtered.push_back(obj);

    }

    return filtered;
        
}

vector<RecJetFormat> filtered_jets(vector<RecJetFormat> jets, float ptmin, float etamax, const vector<RecLeptonFormat>& leptons, float deltaRmin) {
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
        bool passdeltaR = true;

        for (const auto & lepton : leptons) {
            if (jet.dr(lepton) <= deltaRmin) {
                passdeltaR = false;
                break;
            }
        }

        // putting successful jets in the filtered list
        if (passdeltaR) {
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
             std::vector<RecJetFormat>& Ak4jets) {
            
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

    // Sorting jets
    std::sort(Ak4jets.begin(), Ak4jets.end(), leadingjet);

}

// Veto b-tagged jets

bool nobtag(vector<RecJetFormat> jets) {
    for (auto jet : jets) {
        if (jet.btag()) {
            return false;
        }
    }
    return true;
}

// creating the ak15 jet which models the suep

std::pair<std::vector<fastjet::PseudoJet>, std::vector<std::vector<fastjet::PseudoJet>>> getak15jets(std::vector <RecTrackFormat> tracks, std::vector<RecLeptonFormat> leptons, double pt_cut, double eta_cut, double dz_cut, double d0_cut, double dr_cut)
{
    std::vector<fastjet::PseudoJet> input_particles;

    // if(leptons.empty()){
    //     // std::cerr << "ERROR: getak15jets called with empty leptons vector. Skipping." << std::endl;
    //     return {}; // Return empty jets safely
    // }

    // only constructing this jet with particles which pass analysis-level selections

    for (const auto &track : tracks) {
        if (track.pt() >= pt_cut && std::abs(track.eta()) <= eta_cut &&
            std::abs(track.d0()) < d0_cut && std::abs(track.dz()) < dz_cut && 
            track.dr(leptons.at(0)) >= dr_cut
            )

            // creating the particles (constituents) which makeup our jet

            {
            fastjet::PseudoJet particle(track.px(), track.py(), track.pz(), track.e());
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

double deltar(const RecJetFormat& rec_jet, const fastjet::PseudoJet& pseudo_jet)
{
    double rec_eta = rec_jet.eta();
    double rec_phi = rec_jet.phi();

    double pseudo_eta = pseudo_jet.eta();
    double pseudo_phi = pseudo_jet.phi();

    // wrap arounds
    if (pseudo_phi > M_PI) {
        pseudo_phi -= 2 * M_PI;
    }

    // calculate deta and dphi

    double deta = rec_eta - pseudo_eta;
    double dphi = std::atan2(std::sin(rec_phi - pseudo_phi), std::cos(rec_phi - pseudo_phi));

    // calculate deltar

    return std::sqrt(deta * deta + dphi * dphi);
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

    // make the signal region //

    Manager()->AddRegionSelection("SR");

    // add selections //

    Manager()->AddCut("orthogonality");
    Manager()->AddCut("onelep");
    // one ak4
    Manager()->AddCut("oneak15");
    Manager()->AddCut("metptcut");
    Manager()->AddCut("wptcut");
    Manager()->AddCut("onshell");
    Manager()->AddCut("noBs");
    Manager()->AddCut("dphiwsuep");
    Manager()->AddCut("dphimetsuep");
    Manager()->AddCut("dphilepsuep");
    Manager()->AddCut("mindphimetak4cut");
    Manager()->AddCut("minak15ptcut");
    Manager()->AddCut("lepoverlap");
    Manager()->AddCut("overlap");
    Manager()->AddCut("ak4ak15dR");
    Manager()->AddCut("overlap");
    Manager()->AddCut("ratio");
    // delta phi (MET, JET) > 1.5

    // make histos //

    // for the w //
    Manager()->AddHisto("wmass", 60, 0.0, 140.0);
    Manager()->AddHisto("wpt", 300, 0, 300);
    Manager()->AddHisto("weta", 40,-3.14,3.14);
    Manager()->AddHisto("wphi", 40,-3.14,3.14);

    // for the leading lepton //
    Manager()->AddHisto("lep1pt", 300, 0, 300);
    Manager()->AddHisto("lep1eta", 40,-3.14,3.14);
    Manager()->AddHisto("lep1phi", 40,-3.14,3.14);
    Manager()->AddHisto("reconstructed electrons", 10, 0, 10);
    Manager()->AddHisto("looselep", 10, 0, 10);
    Manager()->AddHisto("tightlep", 10, 0, 10);

    // for the leading muon //
    Manager()->AddHisto("muo1pt", 300, 0, 300);
    Manager()->AddHisto("muo1eta", 40,-3.14,3.14);
    Manager()->AddHisto("muo1phi", 40,-3.14,3.14);
    Manager()->AddHisto("reconstructed muons", 10, 0, 10);
    Manager()->AddHisto("loosemu", 10, 0, 10);
    Manager()->AddHisto("tightmu", 10, 0, 10);

    // for the leading electron //
    Manager()->AddHisto("ele1pt", 300, 0, 300);
    Manager()->AddHisto("ele1eta", 40,-3.14,3.14);
    Manager()->AddHisto("ele1phi", 40,-3.14,3.14);
    Manager()->AddHisto("loose_ele", 10, 0, 10);
    Manager()->AddHisto("tightele", 10, 0, 10);



    // for the ak4jets //
    Manager()->AddHisto("NJets", 10,0.0,10.0);
    Manager()->AddHisto("ak41pt", 300, 0, 300);
    Manager()->AddHisto("ak41eta", 40,-3.14,3.14);
    Manager()->AddHisto("ak41phi", 40,-3.14,3.14);
    Manager()->AddHisto("ak41ntracks", 100,0.0,100.0);

    // for the ak15 jets //
    Manager()->AddHisto("ak151pt", 300, 0, 300);
    Manager()->AddHisto("ak151eta", 40,-3.14,3.14);
    Manager()->AddHisto("ak151phi", 40,-3.14,3.14);
    Manager()->AddHisto("ak151ntracks", 100,0.0,100.0);
    Manager()->AddHisto("ak151mass", 400,0,400);

    // abcd region variables & their limits //
    Manager()->AddHisto("ABCD_A", 10, 10.0, 20.0);
    Manager()->AddHisto("ABCD_B", 10, 20.0, 30.0);
    Manager()->AddHisto("ABCD_C", 100, 30.0, 130.0);
    Manager()->AddHisto("ABCD_D", 10, 10.0, 20.0);
    Manager()->AddHisto("ABCD_E", 10, 20.0, 30.0);
    Manager()->AddHisto("ABCD_F0", 10, 30.0, 40.0);
    Manager()->AddHisto("ABCD_F1", 10, 40.0, 50.0);
    Manager()->AddHisto("ABCD_F2", 10, 50.0, 60.0);
    Manager()->AddHisto("ABCD_F3", 20, 60.0, 80.0);
    Manager()->AddHisto("ABCD_F4", 20, 80.0, 100.0);
    Manager()->AddHisto("ABCD_G", 10, 10.0, 20.0);
    Manager()->AddHisto("ABCD_H", 10, 20.0, 30.0);
    Manager()->AddHisto("ABCD_SR0", 10, 30.0, 40.0);
    Manager()->AddHisto("ABCD_SR1", 10, 40.0, 50.0);
    Manager()->AddHisto("ABCD_SR2", 10, 50.0, 60.0);
    Manager()->AddHisto("ABCD_SR3", 10, 60.0, 80.0);
    Manager()->AddHisto("ABCD_SR4", 20, 80.0, 100.0);

    // for met //
    Manager()->AddHisto("metpt", 300, 0, 300);
    Manager()->AddHisto("metphi", 40, -3.14, 3.14);
    Manager()->AddHisto("meteta", 40, -3.14,3.14);

    // sphericity //
    Manager()->AddHisto("sphericity", 50, 0.0, 1.0);
    Manager()->AddHisto("boosted_sphericity", 50, 0.0, 1.0);

    // everything runs smoothly //
    return true;
}

///////////////////////////////////////////////////////////////
//                          execute                          //
//        function called each time an event is read         //
///////////////////////////////////////////////////////////////

bool user::Execute(SampleFormat& sample, const EventFormat& event)
{
    // std::cout << "in the event loop" << std::endl;
    // std::cout << "Number of electrons in event: " << event.rec()->electrons().size() << std::endl;
    // std::cout << "Number of muons in event: " << event.rec()->muons().size() << std::endl;
    // std::cout << "defining object cuts" << std::endl;



    // defining object cuts

    // for electrons 

    float const loose_electron_minpt = 15;
    float const tight_electron_minpt = 35;
    float const electron_maxeta = 2.5;
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

    // for muons

    float const loose_muon_minpt = 10;
    float const tight_muon_minpt = 30;
    float const muon_maxeta = 2.4;
    float const loose_muon_dz = 0.01;
    float const tight_muon_dz = 0.05;
    float const muon_d0 = 0.02;
    // utils ref: https://github.com/SUEPPhysics/SUEPCoffea_dask/blob/d21ebb8d27d8a9c4f264fd680fc2003752aeab8b/workflows/WH_utils.py#L273-L274
    // how to implement: https://github.com/MadAnalysis/madanalysis5/blob/3900ec9002ee7c6963162e72feea58e6971c61eb/tools/SampleAnalyzer/Interfaces/delphes/delphes_cms.tcl#L544

    // for ak4jets

    float const ak4_minpt = 30;
    float const ak4_maxeta = 2.4;
    float const ak4lep_deltar = 0.4;
    // jet id: https://github.com/SUEPPhysics/SUEPCoffea_dask/blob/d21ebb8d27d8a9c4f264fd680fc2003752aeab8b/workflows/WH_utils.py#L67C73-L67C78

    // adding weights
    // std::cout << "adding weights" << std::endl;

    double weight = 1.;
    if (!Configuration().IsNoEventWeight() && event.mc()!=0) {
        weight = event.mc()->weight();
    }

    Manager()->InitializeForNewEvent(weight);
    if (event.rec() == 0) {return true;}

    //////////////////////////////////////////////
    //          applying the selections         //
    //////////////////////////////////////////////


    //////////////////////////////////////////////
    //           object level selections        //
    //////////////////////////////////////////////

    // make loose leptons
    // std::cout << "making loose leptons" << std::endl;


    // setup for isolation

    std::vector<RecTrackFormat> tracks = event.rec()->tracks();
    std::vector<RecParticleFormat> eflowPhotons = event.rec()->EFlowPhotons();
    std::vector<RecParticleFormat> eflowNeutralHadrons = event.rec()->EFlowNeutralHadrons();

    // electrons
    vector<RecLeptonFormat> loose_electrons = filtered_loose_electrons(event.rec()->electrons(), loose_electron_minpt, electron_maxeta, tracks, eflowPhotons, eflowNeutralHadrons);
    // std::cout << "loose_electrons size: " << loose_electrons.size() << std::endl;
    vector<RecLeptonFormat> loose_posElectrons = filtered_loose_electrons(event.rec()->electrons(), loose_electron_minpt, electron_maxeta,tracks, eflowPhotons, eflowNeutralHadrons, "+");
    vector<RecLeptonFormat> loose_negElectrons = filtered_loose_electrons(event.rec()->electrons(), loose_electron_minpt, electron_maxeta, tracks, eflowPhotons, eflowNeutralHadrons, "-");

    // muons
    vector<RecLeptonFormat> loose_muons = filtered_loose_muons(event.rec()->muons(), loose_muon_minpt, muon_maxeta, loose_muon_dz, muon_d0, tracks, eflowPhotons, eflowNeutralHadrons);
    // std::cout << "loose_muons size: " << loose_muons.size() << std::endl;
    vector<RecLeptonFormat> loose_posMuons = filtered_loose_muons(event.rec()->muons(), loose_muon_minpt, muon_maxeta, loose_muon_dz, muon_d0, tracks, eflowPhotons, eflowNeutralHadrons, "+");
    vector<RecLeptonFormat> loose_negMuons = filtered_loose_muons(event.rec()->muons(), loose_muon_minpt, muon_maxeta, loose_muon_dz, muon_d0, tracks, eflowPhotons, eflowNeutralHadrons, "-");


    // leptons (combined electrons and muons)
    vector<RecLeptonFormat> loose_leptons = loose_electrons;
    loose_leptons.insert(loose_leptons.end(), loose_muons.begin(), loose_muons.end());
    // std::cout << "loose_leptons size: " << loose_leptons.size() << std::endl;
    vector<RecLeptonFormat> loose_posLeptons = loose_posElectrons;
    loose_posLeptons.insert(loose_posLeptons.end(), loose_posMuons.begin(), loose_posMuons.end());
    vector<RecLeptonFormat> loose_negLeptons = loose_negElectrons;
    loose_negLeptons.insert(loose_negLeptons.end(), loose_negMuons.begin(), loose_negMuons.end());

    // jets
    vector<RecJetFormat> Ak4jets = filtered_jets(event.rec()->jets(), ak4_minpt, ak4_maxeta, loose_leptons, ak4lep_deltar);

    // sorting loose objects
    leading(loose_electrons, loose_posElectrons, loose_negElectrons, loose_muons, loose_posMuons, loose_negMuons, loose_leptons, loose_posLeptons, loose_negLeptons, Ak4jets);

    //////////////////////////////////////////////
    //           event level selections         //
    //////////////////////////////////////////////
    // std::cout << "doing event level selections" << std::endl;


    // orthogonality to ggf offline: remove events with no muons or electrons
    // std::cout << "ggf orthogonality" << std::endl;



    bool vetonoleptons = (loose_leptons.size() == 0);

    // orthogonality to zh: remove events with a pair of OSSF leptons
    // std::cout << "zh orthogonality" << std::endl;


    bool SF = (loose_muons.size() == 2 || loose_electrons.size() == 2);
    bool OS = (loose_posLeptons.size() == 1 && loose_negLeptons.size() == 1);
    bool OSSF = OS && SF;
    
    // apply both orthogonality cuts
    
    // std::cout << "defining orthogonality cuts" << std::endl;
    bool orthogonality = (!vetonoleptons) && (!OSSF);
    // std::cout << "applying orthogonality cuts" << std::endl;
    if (not Manager()->ApplyCut(orthogonality, "orthogonality")) return true;
    // std::cout << "after orthogonality cuts" << std::endl;

    // apply tight selection on leptons passing orthogonality cuts
    // std::cout << "making tight leptons" << std::endl;
    // std::cout << "loose_muons size before tight: " << loose_muons.size() << std::endl;
    // std::cout << "loose_electrons size before tight: " << loose_electrons.size() << std::endl;
    // std::cout << "tracks size before tight: " << tracks.size() << std::endl;
    // std::cout << "photons size before tight: " << eflowPhotons.size() << std::endl;
    // std::cout << "neutral hadrons size before tight: " << eflowNeutralHadrons.size() << std::endl;


    // electrons
    vector<RecLeptonFormat> tight_electrons = filtered_tight_electrons(loose_electrons, tight_electron_minpt, tracks, eflowPhotons, eflowNeutralHadrons);
    vector<RecLeptonFormat> tight_posElectrons = filtered_tight_electrons(loose_electrons, tight_electron_minpt,tracks, eflowPhotons, eflowNeutralHadrons, "+");
    vector<RecLeptonFormat> tight_negElectrons = filtered_tight_electrons(loose_electrons, tight_electron_minpt, tracks, eflowPhotons, eflowNeutralHadrons, "-");

    // std::cout << "finished with tight electrons" << std::endl;

    // muons
    vector<RecLeptonFormat> tight_muons = filtered_tight_muons(loose_muons, tight_muon_minpt, tight_muon_dz, tracks, eflowPhotons, eflowNeutralHadrons);
    vector<RecLeptonFormat> tight_posMuons = filtered_tight_muons(loose_muons, tight_muon_minpt, tight_muon_dz, tracks, eflowPhotons, eflowNeutralHadrons, "+");
    vector<RecLeptonFormat> tight_negMuons = filtered_tight_muons(loose_muons, tight_muon_minpt, tight_muon_dz, tracks, eflowPhotons, eflowNeutralHadrons, "-");

    // std::cout << "finished with tight muons" << std::endl;


    // leptons (combined electrons and muons)
    vector<RecLeptonFormat> tight_leptons = tight_electrons;
    tight_leptons.insert(tight_leptons.end(), tight_muons.begin(), tight_muons.end());
    vector<RecLeptonFormat> tight_posLeptons = tight_posElectrons;
    tight_posLeptons.insert(tight_posLeptons.end(), tight_posMuons.begin(), tight_posMuons.end());
    vector<RecLeptonFormat> tight_negLeptons = tight_negElectrons;
    tight_negLeptons.insert(tight_negLeptons.end(), tight_negMuons.begin(), tight_negMuons.end());

    // need exactly one tight lepton

    bool onelep = (tight_leptons.size() == 1);
    if (not Manager()->ApplyCut(onelep, "onelep")) return true;

    // std::cout << "finished with tight leptons" << std::endl;
    // std::cout << "tight Leptons size: " << tight_leptons.size() << std::endl;


    // // check sizes again
    // std::cout << "tight_muons size: " << tight_muons.size() << std::endl;
    // std::cout << "tight_electrons size: " << tight_electrons.size() << std::endl;
    // std::cout << "tracks size after tight: " << tracks.size() << std::endl;
    // std::cout << "photons size after tight: " << eflowPhotons.size() << std::endl;
    // std::cout << "neutral hadrons size after tight: " << eflowNeutralHadrons.size() << std::endl;

    // jets
    vector<RecJetFormat> tight_Ak4jets = filtered_jets(event.rec()->jets(), ak4_minpt, ak4_maxeta, tight_leptons, ak4lep_deltar);

    // std::cout << "finished with tight ak4jets" << std::endl;


    // sorting tight objects
    leading(tight_electrons, tight_posElectrons, tight_negElectrons, tight_muons, tight_posMuons, tight_negMuons, tight_leptons, tight_posLeptons, tight_negLeptons, tight_Ak4jets);

    // std::cout << "finished with sorting tight objects" << std::endl;


    // met selections

    float const minmetpt = 30;
    float const mindphimetak4 = 1.5;

    bool metptcut = (event.rec()->MET().pt() >= minmetpt);
    if (!Manager()->ApplyCut(metptcut, "metptcut")) return true;

    // std::cout << "met cut applied" << std::endl;


    // ak4 jets

    if (tight_Ak4jets.size() < 1) {
        // std::cout << "Event rejected: Less than 1 AK4 jet present." << std::endl;
        return true;
    }

    float dphi_ak4_met = fabs(event.rec()->MET().phi() - tight_Ak4jets.at(0).phi());
    if (dphi_ak4_met > M_PI) {
        dphi_ak4_met = 2 * M_PI - dphi_ak4_met;
        }
    bool mindphimetak4cut = (dphi_ak4_met > mindphimetak4);
    if (not Manager()->ApplyCut(mindphimetak4cut, "mindphimetak4cut")) return true;

    // std::cout << "ak4 cuts applied" << std::endl;

    
    // w

    float const minwpt = 60;
    float const minwmass = 30;
    float const maxwmass = 130;

    // tracks

    float const mintrackpt = 1;
    float const maxtracketa = 2.5;
    float const maxtrackd0 = 0.05;
    float const maxtrackdz = 0.05;
    float const mintracklepdr = 0.4;

    // btag veto -- for now

    bool noBs = nobtag(tight_Ak4jets);
    if (not Manager()->ApplyCut(noBs, "noBs")) return true;

    // std::cout << "btag veto applied" << std::endl;


    // ak-15

    float const minak15pt = 60;

    // ak-15 clustering 

    // ak-15 clustering
    // std::cout << "Tracks size: " << event.rec()->tracks().size() << std::endl;
    // std::cout << "Leptons size: " << tight_leptons.size() << std::endl;
    // std::cout << "Mintrackpt: " << mintrackpt << std::endl;
    // std::cout << "Maxtracketa: " << maxtracketa << std::endl;
    // std::cout << "Maxtrackdz: " << maxtrackdz << std::endl;
    // std::cout << "Maxtrackd0: " << maxtrackd0 << std::endl;
    // std::cout << "Mintracklepdr: " << mintracklepdr << std::endl;
    // std::cout << "Calling getak15jets..." << std::endl;
    std::pair<std::vector<fastjet::PseudoJet>, std::vector<std::vector<fastjet::PseudoJet>>> ak15jetsout;
    if (tight_leptons.empty()) {
        // std::cout << "No tight leptons, skipping AK15 clustering." << std::endl;
        return true;
    }
    else {
         ak15jetsout = getak15jets(event.rec()->tracks(), tight_leptons, mintrackpt, maxtracketa, maxtrackdz, maxtrackd0, mintracklepdr);
    }
    
    // std::cout << "getak15jets finished executing." << std::endl;


    // std::cout << "Size of ak15jetsout.first: " << ak15jetsout.first.size() << std::endl;
    // std::cout << "Size of ak15jetsout.second: " << ak15jetsout.second.size() << std::endl;

    // Ensure we have at least 2 elements before accessing index 1
    if (ak15jetsout.second.at(0).size() < 2) {
        std::cerr << "Error: getak15jets returned a vector with insufficient elements (size = " 
                << ak15jetsout.second.at(0).size() << ")" << std::endl;
        return false;
    }

    // Assign values only if the size is sufficient
    // std::cout << "Assigning ak15jets and ak15jetsconst..." << std::endl;
    std::vector<fastjet::PseudoJet> ak15jets = ak15jetsout.first;
    // std::cout << "ak15jets size: " << ak15jets.size() << std::endl;
    std::vector<std::vector<fastjet::PseudoJet>> ak15jetsconst = ak15jetsout.second;
    // std::cout << "ak15jetsconst size: " << ak15jetsconst.size() << std::endl;

    // std::cout << "Printing ak15jets elements:" << std::endl;
    // for (size_t i = 0; i < ak15jets.size(); i++) {
    //     std::cout << "Jet " << i << " pt: " << ak15jets[i].pt() << std::endl;
    // }

    // std::cout << "Printing ak15jetsconst elements:" << std::endl;
    // for (size_t i = 0; i < ak15jetsconst.size(); i++) {
    //     if (ak15jetsconst[i].empty()) {
    //         std::cerr << "WARNING: ak15jetsconst[" << i << "] is empty!" << std::endl;
    //     } else {
    //         std::cout << "ak15jetsconst[" << i << "] size: " << ak15jetsconst[i].size() << std::endl;
    //     }
    // }

    // // need at least one cluster
    bool oneak15 = (ak15jets.size() > 0);
    if (not Manager()->ApplyCut(oneak15, "oneak15")) return true;
    // std::cout << "Checking ak15jetsconst[0] size: " << ak15jetsconst.at(0).size() << std::endl;


   // Check leptons exist
    if (tight_leptons.empty()) {
        std::cerr << "ERROR: No leptons found! Cannot reconstruct W." << std::endl;
        return false;
    }
    // std::cout << "Lepton px: " << tight_leptons.at(0).px() 
    //         << " py: " << tight_leptons.at(0).py() 
    //         << " pz: " << tight_leptons.at(0).pz() 
    //         << " E: " << tight_leptons.at(0).e() << std::endl;

    // Check MET exists
    // if (!event.rec()) {
    //     std::cerr << "ERROR: event.rec() is NULL!" << std::endl;
    //     return false;
    // }
    // std::cout << "MET px: " << event.rec()->MET().px() 
    //         << " py: " << event.rec()->MET().py() 
    //         << " pz: " << event.rec()->MET().pz() 
    //         << " E: " << event.rec()->MET().e() << std::endl;

    // W reconstruction
    ParticleBaseFormat recoW;
    recoW += tight_leptons.at(0).momentum();
    recoW += event.rec()->MET().momentum();

    // std::cout << "W candidate mass: " << recoW.m() << std::endl;
    bool wptcut = (recoW.pt() > minwpt);
    if (not Manager()->ApplyCut(wptcut, "wptcut")) return true;

    // W mass selection
    // std::cout << "Calling ApplyCut for onshell on event " << "..." << std::endl;

    // W mass selection
    // std::cout << "Calling ApplyCut for onshell..." << std::endl;
    bool onshell = (recoW.m() >= minwmass && recoW.m() <= maxwmass);
    // std::cout << "W mass: " << recoW.m() << std::endl;
    // std::cout << "Applying cut: onshell" << std::endl;
    if (!Manager()->ApplyCut(onshell, "onshell")) {
        // std::cout << "Cut failed: onshell" << std::endl;
        return true;
    }


    // ak15 pt
    // std::cout << "Checking if ak15jets is empty..." << std::endl;
    // std::cout << "ak15jets size: " << ak15jets.size() << std::endl;
    if (!ak15jets.empty()) {
        // std::cout << "Leading ak15 jet pt: " << ak15jets.at(0).pt() << std::endl;
    } else {
        std::cerr << "ERROR: ak15jets is empty, cannot access leading jet!" << std::endl;
    }
    bool minak15ptcut = (ak15jets.at(0).pt() > minak15pt);
    // std::cout << "Leading ak15 jet pt: " << ak15jets.at(0).pt() << std::endl;
    if (not Manager()->ApplyCut(minak15ptcut, "minak15ptcut")) return true;

    // ak4 lepton deltaR
    // std::cout << "Checking deltaR for Ak4jets..." << std::endl;
    // std::cout << "Ak4jets size: " << Ak4jets.size() << std::endl;
    // std::cout << "ak15jets size: " << ak15jets.size() << std::endl;
    bool lepoverlap = false;
    if (!tight_Ak4jets.empty() && !ak15jets.empty()) {
        // std::cout << "Calculating deltaR between Ak4jets[0] and ak15jets[0]..." << std::endl;
        double deltaRvalue = deltar(tight_Ak4jets.at(0), ak15jets.at(0));
        // std::cout << "deltaR calculated: " << deltaRvalue << std::endl;
        lepoverlap = (deltaRvalue < 0.4);
    } 
    if (not Manager()->ApplyCut(lepoverlap, "lepoverlap")) return true;

    // ak4 ak15 overlap
    bool overlap = false;
    if (!tight_Ak4jets.empty() && !ak15jets.empty()) {
        // std::cout << "Checking deltaR for Ak4jets and ak15jets..." << std::endl;
        overlap = (deltar(tight_Ak4jets.at(0), ak15jets.at(0)) < 0.4);
    }
    if (not Manager()->ApplyCut(overlap, "overlap")) return true;

    // w_pt / w-suep < 3

    float const maxratio = 3.0;

    bool ratio = (recoW.pt() / ak15jets.at(0).pt() < maxratio);
    if (not Manager()->ApplyCut(ratio, "ratio")) return true;

    // dphi lepton & suep

    float const mindphilepsuep = 1.5;

    float dphi_lep_suep = fabs(tight_leptons.at(0).phi() - ak15jets.at(0).phi());
    if (dphi_lep_suep > M_PI) {
        dphi_lep_suep = 2 * M_PI - dphi_lep_suep;
        }
    bool dphilepsuep = (dphi_lep_suep > mindphilepsuep);
    if (not Manager()->ApplyCut(dphilepsuep, "dphilepsuep")) return true;

    // dphi met & suep

    float const mindphimetsuep = 1.5;

    float dphi_suep_met = fabs(event.rec()->MET().phi() - ak15jets.at(0).phi());
    if (dphi_suep_met > M_PI) {
        dphi_suep_met = 2 * M_PI - dphi_suep_met;
        }
    bool dphimetsuep = (dphi_suep_met > mindphimetsuep);
    if (not Manager()->ApplyCut(dphimetsuep, "dphimetsuep")) return true;

    // dphi w & suep

    float const mindphiwsuep = 1.5;

    float dphi_w_suep = fabs(recoW.phi() - ak15jets.at(0).phi());
    if (dphi_w_suep > M_PI) {
        dphi_w_suep = 2 * M_PI - dphi_w_suep;
        }
    bool dphiwsuep = (dphi_w_suep > mindphiwsuep);
    if (not Manager()->ApplyCut(dphiwsuep, "dphiwsuep")) return true;

    // sphericity

    float const minsphericity = 0.3;

    // bool passed_sphericity_cut = (sphericityval > minsphericity);
    // if (not Manager()->ApplyCut(passed_sphericity_cut, "boosted_sphericity")) return true;


    // boosted
    auto leading_ak15const = ak15jetsout.second.at(0);
    auto leading_ak15jet = ak15jetsout.first.at(0);
    std::vector<fastjet::PseudoJet> boosted_constituents = boost_to_jet_frame(leading_ak15const, leading_ak15jet);
    double boosted_sphericityval = sphericity(boosted_constituents, 1.0);

    // unboosted

    double sphericityval = sphericity(ak15jetsconst.at(0), 1.0);

    // std::cout << "all cuts successfully applied" << std::endl;

    // fill all histos

    // w

    Manager()->FillHisto("wmass", recoW.m());
    Manager()->FillHisto("wpt", recoW.pt());
    Manager()->FillHisto("weta", recoW.eta());
    Manager()->FillHisto("wphi", recoW.phi());

    // muons or electrons, depending on which lepton is in the event



    if(tight_muons.size() > 0) {
        Manager()->FillHisto("muo1pt", tight_muons.at(0).pt());
        Manager()->FillHisto("muo1eta", tight_muons.at(0).eta());
        Manager()->FillHisto("muo1phi", tight_muons.at(0).phi());
        Manager()->FillHisto("tightmu", tight_muons.size());
    }

    else {
        Manager()->FillHisto("ele1pt", tight_electrons.at(0).pt());
        Manager()->FillHisto("ele1eta", tight_electrons.at(0).eta());
        Manager()->FillHisto("ele1phi", tight_electrons.at(0).phi());
        Manager()->FillHisto("tightele", tight_electrons.size());

    }

    if(loose_muons.size() > 0){
        Manager()->FillHisto("loosemu", loose_muons.size());
    }
    else {
        Manager()->FillHisto("loose_ele", loose_electrons.size());
    }

   // fill lepton histos
   Manager()->FillHisto("lep1pt", tight_leptons.at(0).pt());
   Manager()->FillHisto("lep1eta", tight_leptons.at(0).eta());
   Manager()->FillHisto("lep1phi", tight_leptons.at(0).phi());
   Manager()->FillHisto("looselep", loose_leptons.size());
   Manager()->FillHisto("tightlep", tight_leptons.size());


   
    // ak4

    Manager()->FillHisto("NJets", tight_Ak4jets.size());
    Manager()->FillHisto("ak41pt", tight_Ak4jets.at(0).pt());
    Manager()->FillHisto("ak41eta", tight_Ak4jets.at(0).eta());
    Manager()->FillHisto("ak41phi", tight_Ak4jets.at(0).phi());
    Manager()->FillHisto("ak41ntracks", tight_Ak4jets.at(0).ntracks());

    // ak15

    Manager()->FillHisto("ak151pt", ak15jets.at(0).pt());
    Manager()->FillHisto("ak151eta", ak15jets.at(0).eta());
    double phi_recalc = ak15jets.at(0).phi();
    if (phi_recalc > M_PI){
        phi_recalc -= 2 * M_PI;
    }
    Manager()->FillHisto("ak151phi", phi_recalc);
    Manager()->FillHisto("ak151ntracks", ak15jetsconst.at(0).size());
    Manager()->FillHisto("ak151mass", ak15jets.at(0).m());

    // std::cout << "all histos successfully filled" << std::endl;


    // extended abcd regions

    if ((10 <= (ak15jetsconst.at(0).size()) && (ak15jetsconst.at(0).size()) < 20) && (0.3 <= boosted_sphericityval && boosted_sphericityval < 0.4)) {Manager()->FillHisto("ABCD_A", ak15jetsconst.at(0).size());}
    if ((20 <= (ak15jetsconst.at(0).size()) && (ak15jetsconst.at(0).size()) < 30) &&  (0.3 <= boosted_sphericityval && boosted_sphericityval < 0.4)){Manager()->FillHisto("ABCD_B", ak15jetsconst.at(0).size());}
    if ((30 <= (ak15jetsconst.at(0).size()) && (ak15jetsconst.at(0).size()) < 200) &&  (0.3 <= boosted_sphericityval && boosted_sphericityval < 0.4)){Manager()->FillHisto("ABCD_C", ak15jetsconst.at(0).size());}
    if ((10 <= (ak15jetsconst.at(0).size()) && (ak15jetsconst.at(0).size()) < 20) &&  (0.4 <= boosted_sphericityval && boosted_sphericityval < 0.5)){Manager()->FillHisto("ABCD_D", ak15jetsconst.at(0).size());}

    if ((20 <= (ak15jetsconst.at(0).size()) && (ak15jetsconst.at(0).size()) < 30) &&  (0.4 <= boosted_sphericityval && boosted_sphericityval < 0.5)){Manager()->FillHisto("ABCD_E", ak15jetsconst.at(0).size());}
    if ((30 <= (ak15jetsconst.at(0).size()) && (ak15jetsconst.at(0).size()) < 40) &&  (0.4 <= boosted_sphericityval && boosted_sphericityval < 0.5)){Manager()->FillHisto("ABCD_F0", ak15jetsconst.at(0).size());}
    if ((40 <= (ak15jetsconst.at(0).size()) && (ak15jetsconst.at(0).size()) < 50) &&  (0.4 <= boosted_sphericityval && boosted_sphericityval < 0.5)){Manager()->FillHisto("ABCD_F1", ak15jetsconst.at(0).size());}
    if ((50 <= (ak15jetsconst.at(0).size()) && (ak15jetsconst.at(0).size()) < 60) &&  (0.4 <= boosted_sphericityval && boosted_sphericityval < 0.5)){Manager()->FillHisto("ABCD_F2", ak15jetsconst.at(0).size());}
    if ((60 <= (ak15jetsconst.at(0).size()) && (ak15jetsconst.at(0).size()) < 80) &&  (0.4 <= boosted_sphericityval && boosted_sphericityval < 0.5)){Manager()->FillHisto("ABCD_F3", ak15jetsconst.at(0).size());}
    if ((80 <= (ak15jetsconst.at(0).size()) && (ak15jetsconst.at(0).size()) < 200) &&  (0.4 <= boosted_sphericityval && boosted_sphericityval < 0.5)){Manager()->FillHisto("ABCD_F4", ak15jetsconst.at(0).size());}
    
    if ((10 <= (ak15jetsconst.at(0).size()) && (ak15jetsconst.at(0).size()) < 20) &&  (0.5 <= boosted_sphericityval && boosted_sphericityval < 1.0)){Manager()->FillHisto("ABCD_G", ak15jetsconst.at(0).size());}
    if ((20 <= (ak15jetsconst.at(0).size()) && (ak15jetsconst.at(0).size()) < 30) &&  (0.5 <= boosted_sphericityval && boosted_sphericityval < 1.0)){Manager()->FillHisto("ABCD_H", ak15jetsconst.at(0).size());}

    if ((30 <= (ak15jetsconst.at(0).size()) && (ak15jetsconst.at(0).size()) < 40) &&  (0.5 <= boosted_sphericityval && boosted_sphericityval < 1.0)){Manager()->FillHisto("ABCD_SR0", ak15jetsconst.at(0).size());}
    if ((40 <= (ak15jetsconst.at(0).size()) && (ak15jetsconst.at(0).size()) < 50) &&  (0.5 <= boosted_sphericityval && boosted_sphericityval < 1.0)){Manager()->FillHisto("ABCD_SR1", ak15jetsconst.at(0).size());}
    if ((50 <= (ak15jetsconst.at(0).size()) && (ak15jetsconst.at(0).size()) < 60) &&  (0.5 <= boosted_sphericityval && boosted_sphericityval < 1.0)){Manager()->FillHisto("ABCD_SR2", ak15jetsconst.at(0).size());}
    if ((60 <= (ak15jetsconst.at(0).size()) && (ak15jetsconst.at(0).size()) < 80) &&  (0.5 <= boosted_sphericityval && boosted_sphericityval < 1.0)){Manager()->FillHisto("ABCD_SR3", ak15jetsconst.at(0).size());}
    if ((80 <= (ak15jetsconst.at(0).size()) && (ak15jetsconst.at(0).size()) < 200) &&  (0.5 <= boosted_sphericityval && boosted_sphericityval < 1.0)){Manager()->FillHisto("ABCD_SR4", ak15jetsconst.at(0).size());}

    // std::cout << "all abcd successfully made" << std::endl;

    
    // met

    Manager()->FillHisto("metpt", (event.rec)()->MET().pt());
    Manager()->FillHisto("metphi", event.rec()->MET().phi());
    Manager()->FillHisto("meteta", event.rec()->MET().eta());

    // sphericity
    Manager()->FillHisto("sphericity", sphericityval);
    Manager()->FillHisto("boosted_sphericity", boosted_sphericityval);

    // std::cout << "finished" << std::endl;

    return true;
}

///////////////////////////////////////////////////////////////
//                        finalize                           //
//    function called one time at the end of the analysis    //
///////////////////////////////////////////////////////////////

void user::Finalize(const SampleFormat& summary, const std::vector<SampleFormat>& files){}

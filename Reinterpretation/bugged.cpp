#include "SampleAnalyzer/User/Analyzer/user.h"
using namespace MA5;
using namespace std;
#include <random>
#include <Eigen/Dense>


double iso(vector<RecLeptonFormat> leptons, std::vector <RecTrackFormat> tracks, std::vector <RecTrackFormat> eflowPhotons, std::vector <RecTrackFormat> eflowNeutralHadrons, double iso_minpt, std::string iso = "")
{
    float lep1pt = leptons.at(0).pt();
    float deltaRmax = 0.4;
    float ratiomax = 0.0;

    // to be played with
    if (iso == "WP90") {
        ratiomax = 0.25;
        }
    else if (iso == "WP80") {
        ratiomax = 0.15;
        }
    else if (iso == "pfIso2") {
        ratiomax = 0.25;
        }
    else if (iso == "pfIso5") {
        ratiomax = 0.15;
        }

    for (const auto & lepton : leptons) {
        float totalpt = 0.0;
        for (const auto & track : tracks){
            if ((lepton.dr(track) > deltaRmax) && track.pt() > iso_minpt) {
                totalpt += track.pt();
            }
        }
        for (const auto & photon : eflowPhotons){
            if ((lepton.dr(photon) > deltaRmax) && photon.pt() > iso_minpt) {
                totalpt += photon.pt();
            }
        }
        for (const auto & neutral : eflowNeutralHadrons){
            if ((lepton.dr(neutral) > deltaRmax) && neutral.pt() > iso_minpt) {
                totalpt += neutral.pt();
            }
        }
        float pt_ratio = totalpt / lep1pt;
        if (pt_ratio > ratiomax) {
            return 0;
        }
    }
    return 1;
}

double sphericity(const std::vector<fastjet::PseudoJet>& particles, double r) {
    // std::cout << "Inside sphericity(). Input vector size: " << particles.size() << ", r = " << r << std::endl;
    // check for an empty jet

    if (particles.empty()) {
        std::cerr << "Leading jet empty!" << std::endl;
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
            std::cerr << "Warning: Found particle with zero momentum! Skipping." << std::endl;
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
        std::cerr << "Could not calculate eigenvalues..." << std::endl;
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

vector<RecLeptonFormat> loose_muons(const vector<RecLeptonFormat> objects, const float ptmin, const float etamax, const float dz, const float d0) {
    // Selects loose muons that pass analysis-level cuts

    vector<RecLeptonFormat> filtered;

    for ( auto & obj : objects) {
        
        // Charge Selections
        // if ( charge == '+') {
        //     if (obj.charge() <= 0) {
        //         continue;
        //     }
        // } else if ( charge == '-') {
        //     if (obj.charge() >= 0) {
        //         continue;
        //     }
        // }

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

        std::vector<RecTrackFormat> tracks;
        std::vector<RecTrackFormat> eflowPhotons;
        std::vector<RecTrackFormat> eflowNeutralHadrons;

        if (iso({obj}, tracks, eflowPhotons, eflowNeutralHadrons, 0.1, "pfIso2") == 0) {
            continue;
        }

        filtered.push_back(obj);

    }

    return filtered;
        
}

vector<RecLeptonFormat> tight_muons(vector<RecLeptonFormat> objects, float ptmin, float dz, std::string charge = "") {
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
        std::vector<RecTrackFormat> tracks;
        std::vector<RecTrackFormat> eflowPhotons;
        std::vector<RecTrackFormat> eflowNeutralHadrons;

        if (iso({obj}, tracks, eflowPhotons, eflowNeutralHadrons, 0.1, "pfIso5") == 0) {
            continue;
        }

        filtered.push_back(obj);

    }

    return filtered;
        
}

vector<RecLeptonFormat> loose_electrons(vector<RecLeptonFormat> objects, float ptmin, float etamax, std::string charge = "") {
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

        std::vector<RecTrackFormat> tracks;
        std::vector<RecTrackFormat> eflowPhotons;
        std::vector<RecTrackFormat> eflowNeutralHadrons;

        // wp 90

        if (iso({obj}, tracks, eflowPhotons, eflowNeutralHadrons, ptmin, "WP90") == 0) {
            continue;
        }

        filtered.push_back(obj);

    }

    return filtered;
        
}

vector<RecLeptonFormat> tight_electrons(vector<RecLeptonFormat> objects, float ptmin, std::string charge = "") {
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

        std::vector<RecTrackFormat> tracks;
        std::vector<RecTrackFormat> eflowPhotons;
        std::vector<RecTrackFormat> eflowNeutralHadrons;

        
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
    Manager()->AddCut("metptcut");
    Manager()->AddCut("mindphimetak4cut");
    Manager()->AddCut("noBs");
    Manager()->AddCut("oneak15");
    Manager()->AddCut("onshell");
    Manager()->AddCut("wptcut");
    Manager()->AddCut("minak15ptcut");
    Manager()->AddCut("lepoverlap");
    Manager()->AddCut("overlap");
    Manager()->AddCut("wSUEPptRatio");
    Manager()->AddCut("ak4ak15dR");
    Manager()->AddCut("overlap");
    Manager()->AddCut("ratio");
    Manager()->AddCut("dphilepsuep");
    Manager()->AddCut("dphimetsuep");
    Manager()->AddCut("dphiwsuep");
    // Manager()->AddCut("nconst");
    Manager()->AddCut("sphericity");

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

    // for the leading muon //
    Manager()->AddHisto("muo1pt", 300, 0, 300);
    Manager()->AddHisto("muo1eta", 40,-3.14,3.14);
    Manager()->AddHisto("muo1phi", 40,-3.14,3.14);

    // for the leading electron //
    Manager()->AddHisto("ele1pt", 300, 0, 300);
    Manager()->AddHisto("ele1eta", 40,-3.14,3.14);
    Manager()->AddHisto("ele1phi", 40,-3.14,3.14);

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
    Manager()->AddHisto("ABCD_H", 10, 20.0, 30.0);
    Manager()->AddHisto("ABCD_SR0", 10, 30.0, 40.0);
    Manager()->AddHisto("ABCD_SR1", 10, 40.0, 50.0);
    Manager()->AddHisto("ABCD_SR2", 10, 50.0, 60.0);
    Manager()->AddHisto("ABCD_SR3", 10, 60.0, 80.0);
    Manager()->AddHisto("ABCD_SR4", 20, 80.0, 100.0);

    // for met //
    Manager()->AddHisto("metpt", 300, 0, 300);
    Manager()->AddHisto("metphi", 100, -3.14, 3.14);
    Manager()->AddHisto("meteta", 100, -3.14,3.14);

    // sphericity //
    Manager()->AddHisto("sphericity", 50, 0.0, 1.0);

    // everything runs smoothly //
    return true;
}

///////////////////////////////////////////////////////////////
//                          execute                          //
//        function called each time an event is read         //
///////////////////////////////////////////////////////////////

bool user::Execute(SampleFormat& sample, const EventFormat& event)
{
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

    // electrons
    vector<RecLeptonFormat> loose_electrons = loose_electrons(event.rec()->electrons(), loose_electron_minpt, electron_maxeta);
    vector<RecLeptonFormat> loose_posElectrons = loose_electrons(event.rec()->electrons(), loose_electron_minpt, electron_maxeta, "+");
    vector<RecLeptonFormat> loose_negElectrons = loose_electrons(event.rec()->electrons(), loose_electron_minpt, electron_maxeta, "-");

    // muons
    vector<RecLeptonFormat> loose_muons = loose_muons(event.rec()->muons(), loose_muon_minpt, muon_maxeta, loose_muon_dz, muon_d0);
    vector<RecLeptonFormat> loose_posMuons = loose_muons(event.rec()->muons(), loose_muon_minpt, muon_maxeta, loose_muon_dz, muon_d0, "+");
    vector<RecLeptonFormat> loose_negMuons = loose_muons(event.rec()->muons(), loose_muon_minpt, muon_maxeta, loose_muon_dz, muon_d0, "-");


    // leptons (combined electrons and muons)
    vector<RecLeptonFormat> loose_leptons = loose_electrons;
    loose_leptons.insert(loose_leptons.end(), loose_muons.begin(), loose_muons.end());
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

    // orthogonality to ggf offline: remove events with no muons or electrons

    bool vetonoleptons = (loose_leptons.size() == 0);

    // orthogonality to zh: remove events with a pair of OSSF leptons

    bool SF = (loose_muons.size() == 2 || loose_electrons.size() == 2);
    bool OS = (loose_posLeptons.size() == 1 && loose_negLeptons.size() == 1);
    bool OSSF = OS && SF;
    
    // apply both orthogonality cuts
    bool orthogonality = (!vetonoleptons) && (!OSSF);
    if (not Manager()->ApplyCut(orthogonality, "orthogonality")) return true;

    // apply tight selection on leptons passing orthogonality cuts
    /////// change these params based on how you define tight_... functions ///////////

    // electrons
    vector<RecLeptonFormat> tight_electrons = tight_electrons(loose_electrons, tight_electron_minpt);
    vector<RecLeptonFormat> tight_posElectrons = tight_electrons(loose_electrons, tight_electron_minpt, "+");
    vector<RecLeptonFormat> tight_negElectrons = tight_electrons(loose_electrons, tight_electron_minpt, "-");

    // muons
    vector<RecLeptonFormat> tight_muons = tight_muons(loose_muons, tight_muon_minpt, tight_muon_dz);
    vector<RecLeptonFormat> tight_posMuons = tight_muons(loose_muons, tight_muon_minpt, tight_muon_dz, "+");
    vector<RecLeptonFormat> tight_negMuons = tight_muons(loose_muons, tight_muon_minpt, tight_muon_dz, "-");

    // leptons (combined electrons and muons)
    vector<RecLeptonFormat> tight_leptons = tight_electrons;
    loose_leptons.insert(tight_leptons.end(), tight_muons.begin(), tight_muons.end());
    vector<RecLeptonFormat> tight_posLeptons = tight_posElectrons;
    loose_posLeptons.insert(tight_posLeptons.end(), tight_posMuons.begin(), tight_posMuons.end());
    vector<RecLeptonFormat> tight_negLeptons = tight_negElectrons;
    loose_negLeptons.insert(tight_negLeptons.end(), tight_negMuons.begin(), tight_negMuons.end());

    // jets
    vector<RecJetFormat> tight_Ak4jets = filtered_jets(event.rec()->jets(), ak4_minpt, ak4_maxeta, tight_leptons, ak4lep_deltar);

    // sorting tight objects
    leading(tight_electrons, tight_posElectrons, tight_negElectrons, tight_muons, tight_posMuons, tight_negMuons, tight_leptons, tight_posLeptons, tight_negLeptons, tight_Ak4jets);

    // met selections

    float const minmetpt = 30;
    float const mindphimetak4 = 1.5;

    bool metptcut = (event.rec()->MET().pt() >= minmetpt);
    if (!Manager()->ApplyCut(metptcut, "metptcut")) return true;

    // ak4 jets

    if (tight_Ak4jets.size() < 1) {
        std::cout << "Event rejected: Less than 1 AK4 jet present." << std::endl;
        return true;
    }

    float dphi_ak4_met = fabs(event.rec()->MET().phi() - tight_Ak4jets.at(0).phi());
    if (dphi_ak4_met > M_PI) {
        dphi_ak4_met = 2 * M_PI - dphi_ak4_met;
        }
    bool mindphimetak4cut = (dphi_ak4_met > mindphimetak4);
    if (not Manager()->ApplyCut(mindphimetak4cut, "mindphimetak4cut")) return true;
    
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

    // ak-15

    float const minak15pt = 60;

    // ak-15 clustering 

    // ak-15 clustering
    // std::cout << "Tracks size: " << event.rec()->tracks().size() << std::endl;
    // std::cout << "Leptons size: " << leptons.size() << std::endl;
    // std::cout << "Mintrackpt: " << mintrackpt << std::endl;
    // std::cout << "Maxtracketa: " << maxtracketa << std::endl;
    // std::cout << "Maxtrackdz: " << maxtrackdz << std::endl;
    // std::cout << "Maxtrackd0: " << maxtrackd0 << std::endl;
    // std::cout << "Mintracklepdr: " << mintracklepdr << std::endl;
    // std::cout << "Calling getak15jets..." << std::endl;
    auto ak15jetsout = getak15jets(event.rec()->tracks(), tight_leptons, mintrackpt, maxtracketa, maxtrackdz, maxtrackd0, mintracklepdr);
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
    for (size_t i = 0; i < ak15jets.size(); i++) {
        std::cout << "Jet " << i << " pt: " << ak15jets[i].pt() << std::endl;
    }

    // std::cout << "Printing ak15jetsconst elements:" << std::endl;
    for (size_t i = 0; i < ak15jetsconst.size(); i++) {
        if (ak15jetsconst[i].empty()) {
            std::cerr << "WARNING: ak15jetsconst[" << i << "] is empty!" << std::endl;
        } else {
            std::cout << "ak15jetsconst[" << i << "] size: " << ak15jetsconst[i].size() << std::endl;
        }
    }

    // need at least one cluster
    bool oneak15 = (ak15jets.size() > 0);
    if (not Manager()->ApplyCut(oneak15, "oneak15")) return true;
    // std::cout << "Checking ak15jetsconst[0] size: " << ak15jetsconst.at(0).size() << std::endl;
    double sphericityval = sphericity(ak15jetsconst.at(0), 1.0);
    // std::cout << "Sphericity calculated: " << sphericityval << std::endl;

   // Check leptons exist
    if (tight_leptons.empty()) {
        std::cerr << "ERROR: No leptons found! Cannot reconstruct W." << std::endl;
        return false;
    }
    // std::cout << "Lepton px: " << leptons.at(0).px() 
    //         << " py: " << leptons.at(0).py() 
    //         << " pz: " << leptons.at(0).pz() 
    //         << " E: " << leptons.at(0).e() << std::endl;

    // Check MET exists
    if (!event.rec()) {
        std::cerr << "ERROR: event.rec() is NULL!" << std::endl;
        return false;
    }
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

    if (!Manager()) {
        std::cerr << "ERROR: Manager() is null!" << std::endl;
        return false;
    }

    // W mass selection
    // std::cout << "Calling ApplyCut for onshell..." << std::endl;
    bool onshell = (recoW.m() >= minwmass && recoW.m() <= maxwmass);
    // std::cout << "W mass: " << recoW.m() << std::endl;
    // std::cout << "Applying cut: onshell" << std::endl;
    if (!Manager()->ApplyCut(onshell, "onshell")) {
        std::cout << "Cut failed: onshell" << std::endl;
        return true;
    }


    // ak15 pt
    // std::cout << "Checking if ak15jets is empty..." << std::endl;
    // std::cout << "ak15jets size: " << ak15jets.size() << std::endl;
    if (!ak15jets.empty()) {
        std::cout << "Leading ak15 jet pt: " << ak15jets.at(0).pt() << std::endl;
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
        std::cout << "Calculating deltaR between Ak4jets[0] and ak15jets[0]..." << std::endl;
        double deltaRvalue = deltar(tight_Ak4jets.at(0), ak15jets.at(0));
        std::cout << "deltaR calculated: " << deltaRvalue << std::endl;
        lepoverlap = (deltaRvalue < 0.4);
    } else {
        std::cerr << "WARNING: Skipping deltaR because one or both jet vectors are empty!" << std::endl;
    }
    if (not Manager()->ApplyCut(lepoverlap, "lepoverlap")) return true;

    // ak4 ak15 overlap
    bool overlap = false;
    if (!tight_Ak4jets.empty() && !ak15jets.empty()) {
        std::cout << "Checking deltaR for Ak4jets and ak15jets..." << std::endl;
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

    float dphi_w_suep = fabs(recoW.eta() - ak15jets.at(0).phi());
    if (dphi_w_suep > M_PI) {
        dphi_w_suep = 2 * M_PI - dphi_w_suep;
        }
    bool dphiwsuep = (dphi_w_suep > mindphiwsuep);
    if (not Manager()->ApplyCut(dphiwsuep, "dphiwsuep")) return true;

    // sphericity

    float const minsphericity = 0.3;
    bool sphericity = (sphericityval > minsphericity);
    if (not Manager()->ApplyCut(sphericity, "sphericity")) return true;

    
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
    }

    else {
        Manager()->FillHisto("ele1pt", tight_electrons.at(0).pt());
        Manager()->FillHisto("ele1eta", tight_electrons.at(0).eta());
        Manager()->FillHisto("ele1phi", tight_electrons.at(0).phi());
    }

   // fill lepton histos
   Manager()->FillHisto("lep1pt", tight_leptons.at(0).pt());
   Manager()->FillHisto("lep1eta", tight_leptons.at(0).eta());
   Manager()->FillHisto("lep1phi", tight_leptons.at(0).phi());

   
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

    if ((10 <= (ak15jetsconst.at(0).size()) && (ak15jetsconst.at(0).size()) < 20) && (0.3 <= sphericityval && sphericityval < 0.4)) {Manager()->FillHisto("ABCD_A", ak15jetsconst.at(0).size());}
    if ((20 <= (ak15jetsconst.at(0).size()) && (ak15jetsconst.at(0).size()) < 30) &&  (0.3 <= sphericityval && sphericityval < 0.4)){Manager()->FillHisto("ABCD_B", ak15jetsconst.at(0).size());}
    if ((30 <= (ak15jetsconst.at(0).size()) && (ak15jetsconst.at(0).size()) < 200) &&  (0.3 <= sphericityval && sphericityval < 0.4)){Manager()->FillHisto("ABCD_C", ak15jetsconst.at(0).size());}
    if ((10 <= (ak15jetsconst.at(0).size()) && (ak15jetsconst.at(0).size()) < 20) &&  (0.4 <= sphericityval && sphericityval < 0.5)){Manager()->FillHisto("ABCD_D", ak15jetsconst.at(0).size());}

    if ((20 <= (ak15jetsconst.at(0).size()) && (ak15jetsconst.at(0).size()) < 30) &&  (0.4 <= sphericityval && sphericityval < 0.5)){Manager()->FillHisto("ABCD_E", ak15jetsconst.at(0).size());}
    if ((30 <= (ak15jetsconst.at(0).size()) && (ak15jetsconst.at(0).size()) < 40) &&  (0.4 <= sphericityval && sphericityval < 0.5)){Manager()->FillHisto("ABCD_F0", ak15jetsconst.at(0).size());}
    if ((40 <= (ak15jetsconst.at(0).size()) && (ak15jetsconst.at(0).size()) < 50) &&  (0.4 <= sphericityval && sphericityval < 0.5)){Manager()->FillHisto("ABCD_F1", ak15jetsconst.at(0).size());}
    if ((50 <= (ak15jetsconst.at(0).size()) && (ak15jetsconst.at(0).size()) < 60) &&  (0.4 <= sphericityval && sphericityval < 0.5)){Manager()->FillHisto("ABCD_F2", ak15jetsconst.at(0).size());}
    if ((60 <= (ak15jetsconst.at(0).size()) && (ak15jetsconst.at(0).size()) < 80) &&  (0.4 <= sphericityval && sphericityval < 0.5)){Manager()->FillHisto("ABCD_F3", ak15jetsconst.at(0).size());}
    if ((80 <= (ak15jetsconst.at(0).size()) && (ak15jetsconst.at(0).size()) < 200) &&  (0.4 <= sphericityval && sphericityval < 0.5)){Manager()->FillHisto("ABCD_F4", ak15jetsconst.at(0).size());}
    
    if ((10 <= (ak15jetsconst.at(0).size()) && (ak15jetsconst.at(0).size()) < 20) &&  (0.5 <= sphericityval && sphericityval < 1.0)){Manager()->FillHisto("ABCD_G", ak15jetsconst.at(0).size());}
    if ((20 <= (ak15jetsconst.at(0).size()) && (ak15jetsconst.at(0).size()) < 30) &&  (0.5 <= sphericityval && sphericityval < 1.0)){Manager()->FillHisto("ABCD_H", ak15jetsconst.at(0).size());}

    if ((30 <= (ak15jetsconst.at(0).size()) && (ak15jetsconst.at(0).size()) < 40) &&  (0.5 <= sphericityval && sphericityval < 1.0)){Manager()->FillHisto("ABCD_SR0", ak15jetsconst.at(0).size());}
    if ((40 <= (ak15jetsconst.at(0).size()) && (ak15jetsconst.at(0).size()) < 50) &&  (0.5 <= sphericityval && sphericityval < 1.0)){Manager()->FillHisto("ABCD_SR1", ak15jetsconst.at(0).size());}
    if ((50 <= (ak15jetsconst.at(0).size()) && (ak15jetsconst.at(0).size()) < 60) &&  (0.5 <= sphericityval && sphericityval < 1.0)){Manager()->FillHisto("ABCD_SR2", ak15jetsconst.at(0).size());}
    if ((60 <= (ak15jetsconst.at(0).size()) && (ak15jetsconst.at(0).size()) < 80) &&  (0.5 <= sphericityval && sphericityval < 1.0)){Manager()->FillHisto("ABCD_SR3", ak15jetsconst.at(0).size());}
    if ((80 <= (ak15jetsconst.at(0).size()) && (ak15jetsconst.at(0).size()) < 200) &&  (0.5 <= sphericityval && sphericityval < 1.0)){Manager()->FillHisto("ABCD_SR4", ak15jetsconst.at(0).size());}

    // std::cout << "all abcd successfully made" << std::endl;

    
    // met

    Manager()->FillHisto("metpt", event.rec()->MET().pt());
    Manager()->FillHisto("metphi", event.rec()->MET().phi());
    Manager()->FillHisto("meteta", event.rec()->MET().eta());

    // sphericity

    Manager()->FillHisto("sphericity", sphericityval);

    // std::cout << "finished" << std::endl;

    return true;
}

///////////////////////////////////////////////////////////////
//                        finalize                           //
//    function called one time at the end of the analysis    //
///////////////////////////////////////////////////////////////

void user::Finalize(const SampleFormat& summary, const std::vector<SampleFormat>& files){}

#include "SampleAnalyzer/User/Analyzer/user.h"
using namespace MA5;
using namespace std;
#include <random>
#include <Eigen/Dense>

double sphericity(const std::vector<fastjet::PseudoJet>&particles, double r) {

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

        if (p == 0) continue;
        
        // calculating the weight

        double weight = std::pow(p, r - 2.0);

        // calculating the tensor components

        S_xx += px * px * weight;
        S_xy += px * py * weight;
        S_xz += px * pz * weight;
        S_yy += py * py * weight;
        S_yz += py * pz * weight;
        S_zz += pz * pz + weight;

        // calculating the normalization for the components
        norm += std::pow(p, r);

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

    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> solver(S);
    if (solver.info() =! Eigen::Success) {
        std::cerr << "Could not calculate eigenvalues..." std::endl;
        return 0.0;
    }

    // sort the eigenvalues

    Eigen::Vector3d eigenvals = solver.eigenvalues();
    std::sort(eigenvals.data(), eigenvals.data() + eigenvals.size());

    // grab the two smallest eigenvalues

    eigenval1 = eigenvals[0];
    eigenval2 = eigenvals[1];

    // calculate the sphericity

    double sphericityval = 1.5 * (eigenval1 + eigenval2);

    return sphericityval;
}

vector<RecLeptonFormat> filtered_muons(vector<RecLeptonFormat> objects, float ptmin, float etamax, float dz, float d0, std::string charge = "") {
    // Selects muons that pass analysis-level cuts

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

        // eta
        if (fabs(obj.eta()) > etamax) {
            continue;
        }

        // dz

        if (fabs(obj.dz()) > dz) {
            continue;
        }

        // d0

        if (fabs(object.d0()) > d0) {
            continue;
        }

        filtered.push_back(obj);

    }

    return filtered;
        
}

vector<RecLeptonFormat> filtered_electrons(vector<RecLeptonFormat> objects, float ptmin, float etamax, float dz, float d0, std::string charge = "") {
    // Selects electrons that pass analysis-level cuts

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

        // eta
        if (fabs(obj.eta()) > etamax) {
            continue;
        }

        // dz

        if (fabs(obj.dz()) > dz) {
            continue;
        }

        // d0

        if (fabs(object.d0()) > d0) {
            continue;
        }

        filtered.push_back(obj);

    }

    return filtered;
        
}

vector<RecLeptonFormat> filtered_jets(vector<RecLeptonFormat> jets, float ptmin, float etamax, const vector<RecLeptonFormat>& leptons, float deltaRmin) {
    // Selects jets that pass analysis-level cuts

    vector<RecLeptonFormat> filtered;

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
        bool passdeltaR == true;

        for (const auto & lepton : leptons) {
            if (jet.dr() <= deltaRmin) {
                passdeltaR == false
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
             std::vector<RecLeptonFormat>& Ak4jets) {
            
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
    stdd:sort(Ak4jets.begin(), Ak4jets.end(), leadingjet)

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

std::pair<std::vector<fastjet::PseudoJet>, std::vector<std::vector<fastjet::PseudoJet>>> getak15jets(std::vector <recTrackFormat> tracks, std::vector<RecLeptonFormat> leptons, double pt_cut, double eta_cut, double dz_cut, double d0_cut, double dr_cut)
{
    // only constructing this jet with particles which pass analysis-level selections

    for (const auto &tracks : track) {
        if (
            track.pt() >= pt_cut && 
            std::abs(track.eta()) <= eta_cut &&
            track.dr(leptons.at(0)) >= dr_cut &&
            track.dr(leptons.at(1)) >= dr_cut &&
            )

            // creating the particles (constituents) which makeup our jet

            {
            fastjet:PsuedoJet particle(track.px(), track.py(), track.pz(), track.e());
            input_particles.emplace_bark(particle);
            }
    }

    // filling the r = 1.5 jet with its constituents

    fastjet::JetDefinition jet_def(fastjet::antikt_algorithm, 1.5);
    fastjet::ClusterSequence cluster_seq(input_particles, jet_def);

    // again, sort by pt

    std::vector<fastjet::PseudoJet> inclusive_jets = sorted_by_pt(clust_seq.inclusive_jets(0.0));

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
    double pseudo_phi - pseudo_jet.phi();

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
    Manager()->AddCut("oneak15");
    Manager()->AddCut("onshell");
    Manager()->AddCut("wptcut");
    Manager()->AddCut("noBs");
    Manager()->AddCut("minak15ptcut");
    Manager()->AddCut("lepoverlap");

    // make histos //

    // for the w //
    Manager()->AddHisto("wmass", 60, 0.0, 120.0);
    Manager()->AddHisto("wpt", 500, 0.0, 500.0);
    Manager()->AddHisto("weta", 50,-2.5,2.5);
    Manager()->AddHisto("wphi", 50,-2.5,2.5);

    // for the leading lepton //
    Manager()->AddHisto("lep1pt", 500, 0.0, 500.0);
    Manager()->AddHisto("lep1eta", 50,-2.5,2.5);
    Manager()->AddHisto("lep1phi", 50,-2.5,2.5);

    // for the leading muon //
    Manager()->AddHisto("muo1pt", 500, 0.0, 500.0);
    Manager()->AddHisto("muo1eta", 50,-2.5,2.5);
    Manager()->AddHisto("muo1phi", 50,-2.5,2.5);

    // for the leading electron //
    Manager()->AddHisto("ele1pt", 500, 0.0, 500.0);
    Manager()->AddHisto("ele1eta", 50,-2.5,2.5);
    Manager()->AddHisto("ele1phi", 50,-2.5,2.5);

    // for the ak4jets //
    Manager()->AddHisto("NJets", 10,0.0,10.0);
    Manager()->AddHisto("ak41pt", 500,0.0,500.0);
    Manager()->AddHisto("ak41eta", 100,-5.0,5.0);
    Manager()->AddHisto("ak41phi", 40,-3.14,3.14);
    Manager()->AddHisto("ak41ntracks", 100,0.0,100.0);

    // for the ak15 jets //
    Manager()->AddHisto("ak151pt", 500,0.0,500.0);
    Manager()->AddHisto("ak151eta", 100,-5.0,5.0);
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
    Manager()->AddHisto("metpt" 500, 0.0, 500.0);
    Manager()->AddHisto("metphi" 40, -3.14, 3.14);
    Manager()->AddHisto("meteta" 40, -5.0, 5.0);

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

    float const electron_minpt = 35;
    float const electron_maxeta = 2.5;
    // missing: id, and isolation

    // for muons

    float const muon_minpt = 30;
    float const muon_maxeta = 2.4;
    float const d0 = 0.02;
    float const dz = 0.05;
    // missing: id and isolation

    // for ak4jets

    float const ak4_minpt = 30;
    float const ak4_maxeta = 2.4;
    float const ak4lep_deltar = 0.4;
    // missing: id and b-tagging

    // adding weights

    double weight = 1.;
    if (!Configution().IsNoEventWeight(). && event.mc()!=0) {
        weight = event.mc()->weight;
    }

    Manager()->InitializeForNewEvent(weight);
    if (event.rec() == 0) {return true;}

    //////////////////////////////////////////////
    //          applying the selections         //
    //////////////////////////////////////////////


    //////////////////////////////////////////////
    //           object level selections        //
    //////////////////////////////////////////////

    // electrons
    vector<RecLeptonFormat> electrons = filtered_electrons(event.rec()->electrons(), electron_minpt, electron_maxeta);
    vector<RecLeptonFormat> posElectrons = filtered_electrons(event.rec()->electrons(), electron_minpt, electron_maxeta, "+");
    vector<RecLeptonFormat> negeElectrons = filtered_electrons(event.rec()->electrons(), electron_minpt, electron_maxeta, "-");

    // muons
    vector<RecLeptonFormat> muons = filtered_muons(event.rec()->muons(), muon_minpt, muon_maxeta);
    vector<RecLeptonFormat> posMuons = filtered_muons(event.rec()->muons(), muon_minpt, muon_maxeta, "+");
    vector<RecLeptonFormat> negMuons = filtered_muons(event.rec()->muons(), muon_minpt, muon_maxeta, "-");


    // leptons (combined electrons and muons)
    vector<RecLeptonFormat> leptons = electrons;
    leptons.insert(leptons.end(), muons.begin(), muons.end());

    // jets
    vector<RecLeptonFormat> Ak4jets = filtered_electrons(event.rec()->jets(), ak4_minpt, ak4_maxeta, leptons, ak4lep_deltar);

    // sorting all objects
    leading(electrons, posElectrons, negeElectrons, muons, posMuons, negMuons, leptons, Ak4jets);

    //////////////////////////////////////////////
    //           event level selections         //
    //////////////////////////////////////////////

    // orthogonality to ggf offline: remove events with no muons or electrons

    bool vetonoleptons = (leptons.size() == 0);

    // orthogonality to zh: remove events with a pair of OSSF leptons

    bool SF = (muons.size() == 2 || electrons.size() == 2);
    bool OS = (posLeptons.size() == 1 && negLeptons.size() == 1);
    bool OSSF = OS && SF
    
    // apply both orthogonality cuts
    bool orthogonality = vetonoleptons && OSSF;
    if (not Manager()->ApplyCut(orthogonality), "orthogonality") return true;


    // met selections

    float const minmetpt = 30;
    float const mindphimetak4 = 1.5;
    
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
    // missing from pv and puppi weight

    // btag veto -- for now

    bool noBs = nobtag(Ak4jets);
    if (not Manager()->ApplyCut(noBs, "noBs")) return true;

    // ak-15

    float const minak15pt = 60;

    // ak-15 clustering 

    auto ak15jetsout = getak15jets(event.rec()->tracks(), leptons, mintrackpt, maxtracketa, maxtrackdz, maxtrackd0, mintracklepdr);

    std::vector<fastjet::PseudoJet> ak15jets = ak15jetsout.first
    std::vector<std::vector<fastjet::PseudoJet>> ak15jetsconst = ak15jetsout.second

    // need at least one cluster
    bool oneak15 = (ak15jets.size() > 0);
    if (not Manager()->ApplyCut(oneak15), "oneak15") return true;

    double sphericityval = sphericity(ak15jetsconst.at(0), 1.0);

    // W reco

    RecLeptonFormat = leptons.at(0);

    ParticleBaseFormat recoW;
    recoW += leptons.momentum();

    // w mass selection
    
    bool onshell = (recoW.m() >= minwmass && recoW.m() <= maxwmass);
    if (not Manager() -> ApplyCut(onshell), "onshell") return true;

    // w pt selection
    bool wptcut = (recoW.pt() >= minwpt);
    if (not Manager()->ApplyCut(wptcut, "wptcut")) return true;

    // btag veto -- best i could do for now
    bool noBs = 


    // ak15 pt
    bool minak15ptcut = (ak15jets.at(0).pt() > minak15pt);
    if (not Manager()->ApplyCut(minak15ptcut, "minak15ptcut")) return true;

    // ak4 lepton deltar
    bool lepoverlap = false;
    if (!Ak4jets.empty()) {
        droverlap = (deltar(Ak4jets.at(0), leptons.at(0)) < 0.4);
    }
    if (not Manager()->ApplyCut(lepoverlap, "lepoverlap")) return true;

    // ak4 ak15 overlap 
    bool overlap = false;
    if (!Ak4jets.empty()) {
        droverlap = (deltar(Ak4jets.at(0), ak15jets.at(0)) < 0.4);
    }
    if (not Manager()->ApplyCut(overlap, "overlap")) return true;

    // fill all histos

    // w

    Manager()->FillHisto("wmass", recoW.m());
    Manager()->FillHisto("wpt", recoW.pt());
    Manager()->FillHisto("weta", recoW.eta());
    Manager()->FillHisto("wphi", recoW.phi());

    // muons or electrons, depending on which lepton is in the event

    if(muons.size() > 0) {
        Manager()->FillHisto("muo1pt", muons.at(0).pt());
        Manager()->FillHisto("muo1eta", muons.at(0).eta());
        Manager()->FillHisto("muo1phi", muons.at(0).phi());
    }

    else {
        Manager()->FillHisto("ele1pt", electrons.at(0).pt());
        Manager()->FillHisto("ele1eta", electrons.at(0).eta());
        Manager()->FillHisto("ele1phi", electrons.at(0).phi());
    }

    // ak4

    Manager()->AddHisto("NJets", Ak4jets.size());
    Manager()->AddHisto("ak41pt", Ak4jets.at(0).pt());
    Manager()->AddHisto("ak41eta", Ak4jets.at(0).eta());
    Manager()->AddHisto("ak41phi", Ak4jets.at(0).phi());
    Manager()->AddHisto("ak41ntracks", Ak4jets.at(0).ntracks());

    // ak15

    Manager()->FillHisto("ak151pt", Ak15jets.at(0).pt());
    Manager()->FillHisto("ak151eta", Ak15jets.at(0).eta());
    double phi_recalc = Ak15jets.at(0).phi();
    if (phi_recalc > M_PI){
        phi_recalc -= 2 * M_PI
    }
    Manager()->FillHisto("ak151phi", phi_recalc);
    Manager()->FillHisto("ak151ntracks", ak15jetsconst.at(0).size());
    Manager()->FillHisto("ak151mass", Ak15jets.at(0).m());

    // extended abcd regions


    if ((10 <= (ak15jetsconst.at(0).size()) < 20) &&  (0.3 <= sphericityval < 0.4)){Manager()->FillHisto("ABCD_A", ak15jetsconst.at(0).size());}
    if ((20 <= (ak15jetsconst.at(0).size()) < 30) &&  (0.3 <= sphericityval < 0.4)){Manager()->FillHisto("ABCD_B", ak15jetsconst.at(0).size());}
    if ((30 <= (ak15jetsconst.at(0).size()) < 100) &&  (0.3 <= sphericityval < 0.4)){Manager()->FillHisto("ABCD_C", ak15jetsconst.at(0).size());}
    if ((10 <= (ak15jetsconst.at(0).size()) < 20) &&  (0.4 <= sphericityval < 0.5)){Manager()->FillHisto("ABCD_D", ak15jetsconst.at(0).size());}

    if ((20 <= (ak15jetsconst.at(0).size()) < 30) &&  (0.4 <= sphericityval < 0.5)){Manager()->FillHisto("ABCD_E", ak15jetsconst.at(0).size());}
    if ((30 <= (ak15jetsconst.at(0).size()) < 40) &&  (0.4 <= sphericityval < 0.5)){Manager()->FillHisto("ABCD_F0", ak15jetsconst.at(0).size());}
    if ((40 <= (ak15jetsconst.at(0).size()) < 50) &&  (0.4 <= sphericityval < 0.5)){Manager()->FillHisto("ABCD_F1", ak15jetsconst.at(0).size());}
    if ((50 <= (ak15jetsconst.at(0).size()) < 60) &&  (0.4 <= sphericityval < 0.5)){Manager()->FillHisto("ABCD_F2", ak15jetsconst.at(0).size());}
    if ((60 <= (ak15jetsconst.at(0).size()) < 80) &&  (0.4 <= sphericityval < 0.5)){Manager()->FillHisto("ABCD_F3", ak15jetsconst.at(0).size());}
    if ((80 <= (ak15jetsconst.at(0).size()) < 100) &&  (0.4 <= sphericityval < 0.5)){Manager()->FillHisto("ABCD_F4", ak15jetsconst.at(0).size());}
    
    if ((10 <= (ak15jetsconst.at(0).size()) < 20) &&  (0.5 <= sphericityval < 1.0)){Manager()->FillHisto("ABCD_G", ak15jetsconst.at(0).size());}
    if ((20 <= (ak15jetsconst.at(0).size()) < 30) &&  (0.5 <= sphericityval < 1.0)){Manager()->FillHisto("ABCD_H", ak15jetsconst.at(0).size());}

    if ((30 <= (ak15jetsconst.at(0).size()) < 40) &&  (0.5 <= sphericityval < 1.0)){Manager()->FillHisto("ABCD_SR0", ak15jetsconst.at(0).size());}
    if ((40 <= (ak15jetsconst.at(0).size()) < 50) &&  (0.5 <= sphericityval < 1.0)){Manager()->FillHisto("ABCD_SR1", ak15jetsconst.at(0).size());}
    if ((50 <= (ak15jetsconst.at(0).size()) < 60) &&  (0.5 <= sphericityval < 1.0)){Manager()->FillHisto("ABCD_SR2", ak15jetsconst.at(0).size());}
    if ((60 <= (ak15jetsconst.at(0).size()) < 80) &&  (0.5 <= sphericityval < 1.0)){Manager()->FillHisto("ABCD_SR3", ak15jetsconst.at(0).size());}
    if ((80 <= (ak15jetsconst.at(0).size()) < 100) &&  (0.5 <= sphericityval < 1.0)){Manager()->FillHisto("ABCD_SR4", ak15jetsconst.at(0).size());}
    
    // met

    Manager()->FillHisto("metpt" event.rec()->MET().pt());
    Manager()->FillHisto("metphi" event.rec()->MET().phi());
    Manager()->FillHisto("meteta" event.rec()->MET().eta());

    // sphericity

    Manager()->FillHisto("sphericity", sphericityval);

    return true
}

///////////////////////////////////////////////////////////////
//                        finalize                           //
//    function called one time at the end of the analysis    //
///////////////////////////////////////////////////////////////

void user::Finalize(const SampleFormat& summary, const std::vector<SampleFormat>& files)

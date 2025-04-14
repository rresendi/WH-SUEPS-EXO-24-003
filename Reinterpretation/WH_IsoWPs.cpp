#include "SampleAnalyzer/User/Analyzer/user.h"
using namespace MA5;
using namespace std;

bool iso(const RecLeptonFormat &lepton,
         const std::vector<RecTrackFormat> &eflowTracks,
         const std::vector<RecParticleFormat> &eflowPhotons,
         const std::vector<RecParticleFormat> &eflowNeutralHadrons,
         double iso_minpt,
         double deltaRmax,
         double deltaRmin,
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
    for (const auto &track : eflowTracks)
    {
        double dr = lepton.dr(track);
        double pt = track.pt();

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
                                     float etamin,
                                     float etamax,
                                     float iso_pTMin,
                                     float iso_dRMax,
                                     float iso_dRMin,
                                     bool noIso,
                                     const vector<RecTrackFormat> &eflowTracks,
                                     const vector<RecParticleFormat> &eflowPhotons,
                                     const vector<RecParticleFormat> &eflowNeutralHadrons,
                                     const string &selection)
{
    // Helper function to select electrons
    vector<RecLeptonFormat> filtered;

    for (const auto &obj : objects)
    {

        if (obj.pt() < ptmin)
            continue;

        if (fabs(obj.eta()) > etamax)
            continue;

        if (fabs(obj.eta()) < etamin)
            continue;

        if (noIso)
        {
            filtered.push_back(obj);
            continue;
        }
        if (!iso(obj, eflowTracks, eflowPhotons, eflowNeutralHadrons, iso_pTMin, iso_dRMax, iso_dRMin, selection))
        {
            continue;
        }

        filtered.push_back(obj);
    }
    return filtered;
}

vector<RecLeptonFormat> filter_electrons(const vector<RecLeptonFormat> &objects,
                                         float ptmin,
                                         float etamin,
                                         float etamax,
                                         float iso_pTMin,
                                         float iso_dRMax,
                                         float iso_dRMin,
                                         bool noIso,
                                         const vector<RecTrackFormat> &eflowTracks,
                                         const vector<RecParticleFormat> &eflowPhotons,
                                         const vector<RecParticleFormat> &eflowNeutralHadrons,
                                         const string &selection)
{
    // Helper function to select electrons
    vector<RecLeptonFormat> filtered;

    for (const auto &obj : objects)
    {

        if (obj.pt() < ptmin)
            continue;

        if (fabs(obj.eta()) > etamax)
            continue;

        if (fabs(obj.eta()) < etamin)
            continue;

        // Exclude crack region: electrons with 1.444 < |eta| < 1.566
        if ((fabs(obj.eta()) > 1.444) && (fabs(obj.eta()) < 1.566))
            continue;

        if (noIso)
        {
            filtered.push_back(obj);
            continue;
        }
        if (!iso(obj, eflowTracks, eflowPhotons, eflowNeutralHadrons, iso_pTMin, iso_dRMax, iso_dRMin, selection))
        {
            continue;
        }

        filtered.push_back(obj);
    }
    return filtered;
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
    Manager()->AddCut("dummy");

    // ===== Histograms ===== //

    // Muon Histograms
    Manager()->AddHisto("Mu", 50, 0, 50);
    Manager()->AddHisto("looseMu", 50, 0, 50);
    Manager()->AddHisto("vTightMu", 50, 0, 50);

    // Electron Histograms
    Manager()->AddHisto("Ele", 50, 0, 50);
    Manager()->AddHisto("WP90Ele", 50, 0, 50);
    Manager()->AddHisto("WP80Ele", 50, 0, 50);

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
    std::vector<RecTrackFormat> eflowTracks = event.rec()->EFlowTracks();
    std::vector<RecParticleFormat> eflowPhotons = event.rec()->EFlowPhotons();
    std::vector<RecParticleFormat> eflowNeutralHadrons = event.rec()->EFlowNeutralHadrons();

    // ELECTRONS https://twiki.cern.ch/twiki/bin/view/CMS/MultivariateElectronIdentificationRun2
    // Look at the different eta and pt definitions for the iso. Get most of them to be ~consistent with WP90 and WP80 in each eta bin.
    // Just worry about pT > 10
    float const ELECTRON_PT_MIN = 10;
    float const ELECTRON_ETA_MIN = 0;
    float const ELECTRON_ETA_MAX = 0.8;
    float const ELECTRON_ISO_PT_MIN = 0.5;
    float const ELECTRON_ISO_DR_MAX = 0.3;
    float const ELECTRON_ISO_DR_MIN = 0.01;

    // MUONS https://twiki.cern.ch/twiki/bin/viewauth/CMS/SWGuideMuonSelection#Particle_Flow_isolation
    // Mystery???? Try to just use the same pt/eta of the electrons
    float const MUON_PT_MIN = 10;
    float const MUON_ETA_MIN = 0.0;
    float const MUON_ETA_MAX = 0.8;
    float const MUON_ISO_PT_MIN = 0.1;
    float const MUON_ISO_DR_MAX = 0.4;
    float const MUON_ISO_DR_MIN = 0.01;

    //////////////////////////////////////////////
    //  Applying Base Lepton Object Selections  //
    //////////////////////////////////////////////

    // Electron Collections
    vector<RecLeptonFormat> electrons = filter_electrons(event.rec()->electrons(), ELECTRON_PT_MIN, ELECTRON_ETA_MIN, ELECTRON_ETA_MAX, ELECTRON_ISO_PT_MIN, ELECTRON_ISO_DR_MAX, ELECTRON_ISO_DR_MIN, true, eflowTracks, eflowPhotons, eflowNeutralHadrons, "WP90");

    vector<RecLeptonFormat> WP90_electrons = filter_electrons(event.rec()->electrons(), ELECTRON_PT_MIN, ELECTRON_ETA_MIN, ELECTRON_ETA_MAX, ELECTRON_ISO_PT_MIN, ELECTRON_ISO_DR_MAX, ELECTRON_ISO_DR_MIN, false, eflowTracks, eflowPhotons, eflowNeutralHadrons, "WP90");

    vector<RecLeptonFormat> WP80_electrons = filter_electrons(event.rec()->electrons(), ELECTRON_PT_MIN, ELECTRON_ETA_MIN, ELECTRON_ETA_MAX, ELECTRON_ISO_PT_MIN, ELECTRON_ISO_DR_MAX, ELECTRON_ISO_DR_MIN, false, eflowTracks, eflowPhotons, eflowNeutralHadrons, "WP80");

    // Muon Collections

    vector<RecLeptonFormat> muons = filter_muons(event.rec()->muons(), MUON_PT_MIN, MUON_ETA_MIN, MUON_ETA_MAX, MUON_ISO_PT_MIN, MUON_ISO_DR_MAX, MUON_ISO_DR_MIN, true, eflowTracks, eflowPhotons, eflowNeutralHadrons, "pfIso2");

    vector<RecLeptonFormat> looseWP_muons = filter_muons(event.rec()->muons(), MUON_PT_MIN, MUON_ETA_MIN, MUON_ETA_MAX, MUON_ISO_PT_MIN, MUON_ISO_DR_MAX, MUON_ISO_DR_MIN, false, eflowTracks, eflowPhotons, eflowNeutralHadrons, "pfIso2");

    vector<RecLeptonFormat> vTightWP_muons = filter_muons(event.rec()->muons(), MUON_PT_MIN, MUON_ETA_MIN, MUON_ETA_MAX, MUON_ISO_PT_MIN, MUON_ISO_DR_MAX, MUON_ISO_DR_MIN, false, eflowTracks, eflowPhotons, eflowNeutralHadrons, "pfIso5");

    if (not Manager()->ApplyCut(true, "dummy"))
        return true;

    // Muon and Electron Histograms
    Manager()->FillHisto("Mu", muons.size());
    Manager()->FillHisto("looseMu", looseWP_muons.size());
    Manager()->FillHisto("vTightMu", vTightWP_muons.size());

    Manager()->FillHisto("Ele", electrons.size());
    Manager()->FillHisto("WP90Ele", WP90_electrons.size());
    Manager()->FillHisto("WP80Ele", WP80_electrons.size());

    return true;
}

///////////////////////////////////////////////////////////////
//                        Finalize                           //
//    Function called one time at the end of the analysis    //
///////////////////////////////////////////////////////////////

void user::Finalize(const SampleFormat &summary, const std::vector<SampleFormat> &files) {}

double iso(const vector<RecLeptonFormat>& leptons, std::vector <RecTrackFormat> tracks, std::vector <RecTrackFormat> eflowPhotons, std::vector <RecTrackFormat> eflowNeutralHadrons, double minpt, std::string iso = "")
{
    float lep1pt = leptons.at(0).pt();
    bool passdeltaR = true;
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
            if ((lepton.dr.at(0).(track) > deltaRmax) && track.pt() > minpt) {
                totalpt += track.pt();
            }
        }
        for (const auto & photon : eflowPhotons){
            if ((lepton.dr.at(0).(photon) > deltaRmax) && photon.pt() > minpt) {
                totalpt += photon.pt();
            }
        }
        for (const auto & neutral : eflowNeutralHadrons){
            if ((lepton.dr.at(0).(neutral) > deltaRmax) && neutral.pt() > minpt) {
                totalpt += neutral.pt();
            }
        }

    }

    float pt_ratio = totalpt/lep1pt;
    if (pt_ratio > ratiomax){
        return 0;
    }
    return 1;
}

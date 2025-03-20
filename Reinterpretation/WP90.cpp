double WP90(const vector<RecLeptonFormat>& leptons, std::vector <RecTrackFormat> tracks, std::vector <RecTrackFormat> eflowPhotons, std::vector <RecTrackFormat> eflowNeutralHadrons, double minpt)
{
    float lep1pt = leptons.at(0).pt();
    bool passdeltaR = true;
    float deltaRmax = 0.4;
    float ratiomax = 0.25; // to be played with

    for (const auto & lepton : leptons) {
        float totalpt = 0.0;
        for (const auto & track : tracks){
            if ((lepton.at(0).dr(track) > deltaRmax) && track.pt() > minpt) {
                totalpt += track.pt();
            }
        }
        for (const auto & photon : eflowPhotons){
            if ((lepton.at(0).dr(photon) > deltaRmax) && photon.pt() > minpt) {
                totalpt += photon.pt();
            }
        }
        for (const auto & neutral : eflowNeutralHadrons){
            if ((lepton.at(0).dr(neutral) > deltaRmax) && neutral.pt() > minpt) {
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

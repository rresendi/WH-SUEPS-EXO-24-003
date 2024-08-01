def load_lumi_mask(file_path):
    with open(file_path, 'r') as f:
        goldenJSONDict = json.load(f)

    def mask(run, luminosityBlock):
        mask_array = np.zeros_like(run, dtype=bool)
        for i in range(len(run)):
            run_str = str(run[i])
            if run_str in goldenJSONDict:
                goodIntervals = goldenJSONDict[run_str]
                for interval in goodIntervals:
                    min_lumi, max_lumi = interval
                    if min_lumi <= luminosityBlock[i] <= max_lumi:
                        mask_array[i] = True
                        break
        return mask_array

    return mask

def reference_cuts(events, era):
    # Apply lumi mask
    if era == "2016" or era == "2016apv":
        LumiJSON = load_lumi_mask("eos/user/r/rresendi/Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON.txt")
    elif era == "2016apv":
        LumiJSON = load_lumi_mask("eos/user/r/rresendi/Cert_271036-284044_13TeV_Legacy2016APV_JSON_scout.txt")
    elif era == "2017":
        LumiJSON = load_lumi_mask("eos/user/r/rresendi/Cert_294927-306462_13TeV_UL2017_Collisions17_GoldenJSON.txt")
    elif era == "2018":
        LumiJSON = load_lumi_mask("eos/user/r/rresendi/Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.txt")
    else:
        raise ValueError("No era is defined. Please specify the year")

    lumi_mask = LumiJSON(events["run"], events["luminosityBlock"])
    events = events[lumi_mask]

    # Apply the MET cut
    met_cut = events["MET_pt"] > 150
    events = events[met_cut]

    # Apply quality filters
    if era == "2018" or era == "2017":
        cutAnyFilter = (
            (events["Flag_goodVertices"])
            & (events["Flag_globalSuperTightHalo2016Filter"])
            & (events["Flag_HBHENoiseFilter"])
            & (events["Flag_HBHENoiseIsoFilter"])
            & (events["Flag_EcalDeadCellTriggerPrimitiveFilter"])
            & (events["Flag_BadPFMuonFilter"])
            & (events["Flag_BadPFMuonDzFilter"])
            & (events["Flag_eeBadScFilter"])
            & (events["Flag_ecalBadCalibFilter"])
        )
    elif era == "2016" or era == "2016apv":
        cutAnyFilter = (
            (events["Flag_goodVertices"])
            & (events["Flag_globalSuperTightHalo2016Filter"])
            & (events["Flag_HBHENoiseFilter"])
            & (events["Flag_HBHENoiseIsoFilter"])
            & (events["Flag_EcalDeadCellTriggerPrimitiveFilter"])
            & (events["Flag_BadPFMuonFilter"])
            & (events["Flag_BadPFMuonDzFilter"])
            & (events["Flag_eeBadScFilter"])
        )
    else:
        raise ValueError("Unsupported era for quality filters")

    events = events[cutAnyFilter]

    return events

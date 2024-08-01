# grab events that pass reference cuts
def reference_cuts(events, self, met):
    # Apply lumi mask
    if (self.era == "2016" or self.era == "2016apv"):
        LumiJSON = lumi_tools.LumiMask(
            "data/GoldenJSON/Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON.txt"
        )
    elif self.era == "2016":
        LumiJSON = lumi_tools.LumiMask(
            "data/GoldenJSON/Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON_scout.txt"
        )
    elif self.era == "2016apv":
        LumiJSON = lumi_tools.LumiMask(
            "data/GoldenJSON/Cert_271036-284044_13TeV_Legacy2016_Collisions16APV_JSON_scout.txt"
        )
    elif self.era == "2017":
        LumiJSON = lumi_tools.LumiMask(
            "data/GoldenJSON/Cert_294927-306462_13TeV_UL2017_Collisions17_GoldenJSON.txt"
        )
    elif self.era == "2018":
        LumiJSON = lumi_tools.LumiMask(
            "data/GoldenJSON/Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.txt"
        )
    else:
        raise ValueError("No era is defined. Please specify the year")
        
    events = events[LumiJSON(events.run, events.luminosityBlock)]

    # Apply the MET cut
    met_cut = events.MET.pt > 150
    events = events[met_cut]

    # Apply quality filters
    if self.era == "2018" or self.era == "2017":
        cutAnyFilter = (
            (events.Flag.goodVertices)
            & (events.Flag.globalSuperTightHalo2016Filter)
            & (events.Flag.HBHENoiseFilter)
            & (events.Flag.HBHENoiseIsoFilter)
            & (events.Flag.EcalDeadCellTriggerPrimitiveFilter)
            & (events.Flag.BadPFMuonFilter)
            & (events.Flag.BadPFMuonDzFilter)
            & (events.Flag.eeBadScFilter)
            & (events.Flag.ecalBadCalibFilter)
        )
    elif self.era == "2016" or self.era == "2016apv":
        cutAnyFilter = (
            (events.Flag.goodVertices)
            & (events.Flag.globalSuperTightHalo2016Filter)
            & (events.Flag.HBHENoiseFilter)
            & (events.Flag.HBHENoiseIsoFilter)
            & (events.Flag.EcalDeadCellTriggerPrimitiveFilter)
            & (events.Flag.BadPFMuonFilter)
            & (events.Flag.BadPFMuonDzFilter)
            & (events.Flag.eeBadScFilter)
        )
    else:
        raise ValueError("Unsupported era for quality filters")

    events = events[cutAnyFilter]

    return events

import uproot
import os
import argparse
import awkward as ak
import numpy as np
import ROOT
from array import array

# Sets batch mode so no popup window                                                                                                                                                                       
ROOT.gROOT.SetBatch(True)

# Initialize parser                                                                                                                                                                                        
parser = argparse.ArgumentParser()

parser.add_argument("--input", help = "Name of input file", type = str)
args = vars(parser.parse_args())

# Name of sample                                                                                                                                                                                           
sample_name = args["input"]

output_file = "MC_electron_efficiencies.root"
input_file = "/eos/user/j/jreicher/SUEP/WH_private_signals/merged/" + sample_name + ".root"

# suep decay type                                                                                                                                                                                          

if "generic" in sample_name:
    decay_type = "generic"
elif "hadronic" in sample_name:
    decay_type = "hadronic"
else:
    decay_type = "leptonic"

# conditions for what year                                                                                                                                                                                 

if "UL18" in sample_name:
    year = "2018 conditions"
    folder = "ele_eff_outputs_2018/"
elif "UL17" in sample_name:
    year = "2017 conditions"
    folder = "ele_eff_outputs_2017/"
elif "UL16APV" in sample_name:
    folder = "ele_eff_outputs_2016APV/"
else:
    year = "2016 conditions"
    folder = "ele_eff_outputs_2016/"

# dark meson (phi) mass                                                                                                                                                                                    

if "MD2.00" in sample_name:
    md = "2.00 [GeV]"
elif "MD4.00" in sample_name:
    md = "4.00 [GeV]"
elif "MD3.00" in sample_name:
    md = "3.00 [GeV]"
elif "MD8.00" in sample_name:
    md = "8.00 [GeV]"
elif "MD1.00" in sample_name:
    md = "1.00 [GeV]"
else:
    md = "1.40 [GeV]"

# temperature                                                                                                                                                                                              

if "T0.25" in sample_name:
    temp = "0.25"
if "T0.35" in sample_name:
    temp = "0.35"
if "T0.50" in sample_name:
    temp = "0.50"
elif "T0.75" in sample_name:
    temp = "0.75"
elif "T1.00" in sample_name:
    temp = "1.00"
elif "T1.50" in sample_name:
    temp = "1.50"
elif "T2.00" in sample_name:
    temp = "2.00"
elif "T3.00" in sample_name:
    temp = "3.00"
elif "T4.00" in sample_name:
    temp = "4.00"
elif "T8.00" in sample_name:
    temp = "8.00"
elif "T12.00" in sample_name:
    temp = "12.00"
elif "T16.00" in sample_name:
    temp = "16.00"
elif "T32.00" in sample_name:
    temp = "32.00"
else:
    temp = "6.00"

# Gets relevant events from file                                                                                                                                                                           

def Events(f):
    evs = f['Events'].arrays(['HLT_Ele32_WPTight_Gsf',
                'HLT_Ele115_CaloIdVT_GsfTrkIdT',
                'HLT_Photon175',
                'HLT_Photon200',
                'Electron_cutBased',
                'Electron_pt',
                'Electron_mvaFall17V2Iso_WP80',
                'Electron_eta',
                'Electron_dxy',
                'Electron_dz',
                'Electron_phi',
                'Electron_mass',
                'Electron_pdgId',
                'TrigObj_pt',
                'TrigObj_eta',
                'TrigObj_phi',
                'TrigObj_id',
                'TrigObj_filterBits'
                ])
    return evs

with uproot.open(input_file) as f:
    evs = Events(f)
    
# Defining an offline electron                                                                                                                                                                                 
electrons = ak.zip({
        "pt": evs["Electron_pt"],
        "eta": evs["Electron_eta"],
        "phi": evs["Electron_phi"],
        "mass": evs["Electron_mass"],
        "charge": evs["Electron_pdgId"]/(-11),
        "pdgId": evs["Electron_pdgId"],
        "isTight": evs["Electron_mvaFall17V2Iso_WP80"],
        "isTightIso": evs["Electron_mvaFall17V2Iso_WP80"]
        }, with_name = "Momentum4D")

cutElectrons = (
        (evs["Electron_cutBased"] >= 2)
        & (evs["Electron_mvaFall17V2Iso_WP80"])
        & (abs(evs["Electron_dxy"]) < 0.05 + 0.05 * (abs(evs["Electron_eta"]) > 1.479))
        & (abs(evs["Electron_dz"]) < 0.10 + 0.10 * (abs(evs["Electron_eta"]) > 1.479))
        & ((abs(evs["Electron_eta"]) < 1.444) | (abs(evs["Electron_eta"]) > 1.566))
        & (abs(evs["Electron_eta"]) < 2.5)
    )

offlineElectrons = electrons[cutElectrons]

# Gets matched online/offline electrons from file
def isHLTMatched(events, offlineElectrons):
    trigObj = ak.zip({
            "pt": events["TrigObj_pt"],
            "eta": events["TrigObj_eta"],
            "phi": events['TrigObj_phi'],
            "mass": 0.,
            "id": events['TrigObj_id'],
            "filterBits": events['TrigObj_filterBits']
            }, with_name = "Momentum4D")
    
    # Defining the conditions for filtering each trigger
    
    filterbits1 = ((events['TrigObj_filterBits'] & 2) == 2)
    filterbits2 = ((events['TrigObj_filterBits'] & 2048) == 2048)
    filterbits3 = ((events['TrigObj_filterBits'] & 8192) == 8192)
    
    trigObjSingleEl = trigObj[((abs(trigObj.id) == 11)
                               & (trigObj.pt >= 35)
                               & ((abs(trigObj.eta) < 1.444) | (abs(trigObj.eta) > 1.566))
                               & (abs(trigObj.eta) < 2.5)
                               & (filterbits1 |
                                filterbits2 |
                                filterbits3))]
    
    # Computes deltaR2                                                                                                                                                                                         

    def deltaR2(eta1, phi1, eta2, phi2):
        deta = eta1 - eta2
        dphi = phi1 - phi2
        dphi = np.arccos(np.cos(dphi))

        return deta**2 + dphi**2
    
    toMatch1El, trigObjSingleEl = ak.unzip(ak.cartesian([offlineElectrons, trigObjSingleEl], axis=1, nested = True))
    alldr2 = deltaR2(toMatch1El.eta, toMatch1El.phi, trigObjSingleEl.eta, trigObjSingleEl.phi)
    min_dr2 = ak.min(alldr2, axis=2)
    match1El = ak.any(min_dr2 < 0.1, axis=1)
    
    return match1El, min_dr2

# Defines binning and histograms
ele_bin_edges = array('d',[0,2,4,6,8,10,12,
                         14,16,18,20,22,
                         24,26,28,30,32,
                         34,36,38,40,50,
                         60,70,80,90,100,
                         120,140,160,180,200])

# Histograms for overall efficiency
ele_totalhist = ROOT.TH1D("total_events","Total Events",len(ele_bin_edges)-1,ele_bin_edges)
ele_filthist = ROOT.TH1D("filt_events","Filtered Events",len(ele_bin_edges)-1,ele_bin_edges)

# Split into three regions of eta
eta1_ele_totalhist = ROOT.TH1D("total_events","Total Events",len(ele_bin_edges)-1,ele_bin_edges)
eta1_ele_filthist = ROOT.TH1D("filt_events","Filtered Events",len(ele_bin_edges)-1,ele_bin_edges)
eta2_ele_totalhist = ROOT.TH1D("total_events","Total Events",len(ele_bin_edges)-1,ele_bin_edges)
eta2_ele_filthist = ROOT.TH1D("filt_events","Filtered Events",len(ele_bin_edges)-1,ele_bin_edges)
eta3_ele_totalhist = ROOT.TH1D("total_events","Total Events",len(ele_bin_edges)-1,ele_bin_edges)
eta3_ele_filthist = ROOT.TH1D("filt_events","Filtered Events",len(ele_bin_edges)-1,ele_bin_edges)

def ele_hists(events, etas, hists):
    ele_totalhist = hists[0]
    ele_filthist = hists[1]
    eta_min = etas[0]
    eta_max = etas[1]

    # Electron selection                                                                                                                                                                                   

    ele_quality_check = isHLTMatched(events, offlineElectrons)

    # Trigger selection                                                                                                                                                                                    

    if "UL17" in sample_name or "UL18" in sample_name:
        triggerSingleElectron = (
            events["HLT_Ele32_WPTight_Gsf"] |
            events["HLT_Ele115_CaloIdVT_GsfTrkIdT"] |
            events["HLT_Photon200"]
        )
    else:
        triggerSingleElectron = (
            events["HLT_Ele32_WPTight_Gsf"] |
            events["HLT_Ele115_CaloIdVT_GsfTrkIdT"] |
            events["HLT_Photon175"]
        )

    # Cut on eta                                                                                                                                                                                           

    eta_split = (
        (np.abs(events["Electron_eta"]) >= eta_min)
        & (np.abs(events["Electron_eta"]) < eta_max)
    )

    # Select based on trigger                                                                                                                                                                              

    ele = events["Electron_pt"]
    evs = ele[ele_quality_check & eta_split]
    tr_evs = evs[triggerSingleElectron]

    # Fill histograms                                                                                                                                                                                      

    for ev in evs:
        for entry in ev:
            ele_totalhist.Fill(entry)
    for ev in tr_evs:
        for entry in ev:
            ele_filthist.Fill(entry)

    return 0

eta_split = [[0, 1.0], [1.0, 2.0], [2.0, 3.0]]
colors = [ROOT.kBlack, ROOT.kRed, ROOT.kBlue]

c1 = ROOT.TCanvas("c1", "Electron pT vs Min DeltaR", 800, 600)
legend = ROOT.TLegend(0.5, 0.1, 0.9, 0.4)
graphs = []

for i, (eta_min, eta_max) in enumerate(eta_split):
    eta_mask = (np.abs(evs["Electron_eta"]) >= eta_min) & (np.abs(evs["Electron_eta"]) < eta_max)
    matched_electrons = offlineElectrons[eta_mask]
    match1El, min_dr2 = isHLTMatched(evs, matched_electrons)

    matched_pts = matched_electrons.pt[match1El].to_numpy().flatten()
    min_deltaRs = np.sqrt(min_dr2[match1El].to_numpy().flatten())

    graph = ROOT.TGraph(len(matched_pts), array('d', matched_pts), array('d', min_deltaRs))
    graph.SetMarkerStyle(20)
    graph.SetMarkerColor(colors[i])
    graph.SetLineColor(colors[i])
    graphs.append(graph)

    legend.AddEntry(graph, f"{eta_min}<|#eta|<{eta_max}", "l")

graphs[0].SetTitle("Electron pT vs Min DeltaR;Electron pT [GeV];Min DeltaR")
graphs[0].Draw("AP")
for graph in graphs[1:]:
    graph.Draw("P same")

legend.Draw()
c1.SaveAs(sample_name + "_pt_vs_minDeltaR.pdf")

print("Sample " + sample_name + " processed")

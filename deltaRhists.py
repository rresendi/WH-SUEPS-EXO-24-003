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
parser.add_argument("--input", help="Name of input file", type=str)
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

# Defining a good electron                                                                                                                                                                                 

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
        & (evs["Electron_pt"] >= 15)
        & (evs["Electron_mvaFall17V2Iso_WP80"])
        & (abs(evs["Electron_dxy"]) < 0.05 + 0.05 * (abs(evs["Electron_eta"]) > 1.479))
        & (abs(evs["Electron_dz"]) < 0.10 + 0.10 * (abs(evs["Electron_eta"]) > 1.479))
        & ((abs(evs["Electron_eta"]) < 1.444) | (abs(evs["Electron_eta"]) > 1.566))
        & (abs(evs["Electron_eta"]) < 2.5)
    )

goodElectrons = electrons[cutElectrons]

# Computes deltaR2
def deltaR2(eta1, phi1, eta2, phi2):
    deta = eta1 - eta2
    dphi = np.arccos(np.cos(phi1 - phi2))

    return deta ** 2 + dphi ** 2

# Gets matched online/offline electrons from file
def isHLTMatched(events, electrons):
    trigObj = ak.zip({
        "pt": events["TrigObj_pt"],
        "eta": events["TrigObj_eta"],
        "phi": events['TrigObj_phi'],
        "mass": 0.,
        "id": events['TrigObj_id'],
        "filterBits": events['TrigObj_filterBits']
    }, with_name="Momentum4D")

    trigObjSingleEl = trigObj[((abs(trigObj.id) == 11) &
                               ((events['TrigObj_filterBits'] & 2) |
                                (events['TrigObj_filterBits'] & 2048) |
                                (events['TrigObj_filterBits'] & 8192)))]

    toMatch1El, trigObjSingleEl = ak.unzip(ak.cartesian([goodElectrons, trigObjSingleEl], axis=1, nested=True))
    alldr2 = deltaR2(toMatch1El.eta, toMatch1El.phi, trigObjSingleEl.eta, trigObjSingleEl.phi)
    match1El = (ak.sum(ak.where(ak.min(alldr2, axis=2) < 0.1, True, False), axis=1) >= 1)

    return match1El

# Defines binning and histograms
ele_bin_edges = array('d', [0, 2, 4, 6, 8, 10, 12,
                            14, 16, 18, 20, 22,
                            24, 26, 28, 30, 32,
                            34, 36, 38, 40, 50,
                            60, 70, 80, 90, 100,
                            120, 140, 160, 180, 200])

# Create a 2D histogram for deltaR vs. pT for each eta bin
deltaR_vs_pt_eta1 = ROOT.TH2F("deltaR_vs_pt_eta1", "DeltaR vs pT (|eta| < 1)", 30, 0, 200, 100, 0, 0.5)
deltaR_vs_pt_eta2 = ROOT.TH2F("deltaR_vs_pt_eta2", "DeltaR vs pT (1 <= |eta| < 2)", 30, 0, 200, 100, 0, 0.5)
deltaR_vs_pt_eta3 = ROOT.TH2F("deltaR_vs_pt_eta3", "DeltaR vs pT (2 <= |eta| < 2.5)", 30, 0, 200, 100, 0, 0.5)

# Function for filling the deltaR histograms and 2D histograms
def ele_hists(events, etas, deltaR_vs_pt_hist):
    eta_min, eta_max = etas

    # Electron selection
    ele_quality_check = isHLTMatched(events, goodElectrons)

    # Cut on eta
    eta_split = (
        (np.abs(events["Electron_eta"]) >= eta_min) &
        (np.abs(events["Electron_eta"]) < eta_max)
    )

    # Compute deltaR and fill 2D histograms
    trigObj = ak.zip({
        "eta": events["TrigObj_eta"],
        "phi": events['TrigObj_phi']
    }, with_name="Momentum4D")
    trigObjSingleEl = trigObj[((abs(events['TrigObj_id']) == 11) &
                               ((events['TrigObj_filterBits'] & 2) |
                                (events['TrigObj_filterBits'] & 2048) |
                                (events['TrigObj_filterBits'] & 8192)))]
    toMatch1El, trigObjSingleEl = ak.unzip(ak.cartesian([goodElectrons, trigObjSingleEl], axis=1, nested=True))
    alldr2 = deltaR2(toMatch1El.eta, toMatch1El.phi, trigObjSingleEl.eta, trigObjSingleEl.phi)
    deltaR_values = ak.min(alldr2, axis=2)
    electron_pts = events["Electron_pt"]

    for deltaR, pt in zip(deltaR_values[ele_quality_check & eta_split], electron_pts[ele_quality_check & eta_split]):
        for value, pt_val in zip(deltaR, pt):
            deltaR_vs_pt_hist.Fill(pt_val, value)

# Fill 2D histograms
eta_bins = [(0, 1, deltaR_vs_pt_eta1), (1, 2, deltaR_vs_pt_eta2), (2, 2.5, deltaR_vs_pt_eta3)]
for (eta_min, eta_max, deltaR_vs_pt_hist) in eta_bins:
    ele_hists(evs, (eta_min, eta_max), deltaR_vs_pt_hist)

# Create canvas for plotting
c1 = ROOT.TCanvas("canvas", "", 800, 600)
c1.SetRightMargin(0.15)  # Adjust right margin to accommodate legend

# Draw the first 2D histogram
deltaR_vs_pt_eta1.Draw("COLZ")
deltaR_vs_pt_eta1.GetXaxis().SetTitle("pT [GeV]")
deltaR_vs_pt_eta1.GetYaxis().SetTitle("DeltaR")
c1.Update()

# Draw subsequent histograms on top of the first one
deltaR_vs_pt_eta2.Draw("COLZ SAME")
deltaR_vs_pt_eta3.Draw("COLZ SAME")

# Add a legend
legend = ROOT.TLegend(0.7, 0.15, 0.85, 0.4)
legend.AddEntry(deltaR_vs_pt_eta1, "|eta| < 1", "l")
legend.AddEntry(deltaR_vs_pt_eta2, "1 <= |eta| < 2", "l")
legend.AddEntry(deltaR_vs_pt_eta3, "2 <= |eta| < 2.5", "l")
legend.SetTextSize(0.03)
legend.Draw()

# Save the plot
c1.SaveAs(folder + sample_name + "_DeltaR_vs_pT_by_Eta.pdf")

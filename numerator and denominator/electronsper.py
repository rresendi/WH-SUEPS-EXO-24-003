import uproot
import argparse
import numpy as np
import ROOT
import awkward as ak
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
decay_type = ""
if "generic" in sample_name:
    decay_type = "generic"
elif "hadronic" in sample_name:
    decay_type = "hadronic"
else:
    decay_type = "leptonic"

# conditions for what year
year = ""
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
md = ""
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
temp = ""
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
                ])
    return evs

with uproot.open(input_file) as f:
    evs = Events(f)

# Defines binning and histograms
ele_bin_edges = array('d', [0, 2, 4, 6, 8, 10, 12,
                            14, 16, 18, 20, 22,
                            24, 26, 28, 30, 32,
                            34, 36, 38, 40, 50,
                            60, 70, 80, 90, 100,
                            120, 140, 160, 180, 200])

# Histograms for sums in eta regions
eta_ele_sumhist1 = ROOT.TH1D("eta1_sum_events", "Sum Events (|eta|<1.0)", len(ele_bin_edges) - 1, ele_bin_edges)
eta_ele_sumhist2 = ROOT.TH1D("eta2_sum_events", "Sum Events (1.0<|eta|<2.0)", len(ele_bin_edges) - 1, ele_bin_edges)
eta_ele_sumhist3 = ROOT.TH1D("eta3_sum_events", "Sum Events (2.0<|eta|<3.0)", len(ele_bin_edges) - 1, ele_bin_edges)

# Function for filling the histograms
def ele_sums(events, etas, hist):
    eta_min = etas[0]
    eta_max = etas[1]

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

    # quality requirements for electrons
    ele_quality_check = (
            (events["Electron_cutBased"] >= 2)
            & (events["Electron_mvaFall17V2Iso_WP80"])
            & (abs(events["Electron_dxy"]) < 0.05 + 0.05 * (abs(events["Electron_eta"]) > 1.479))
            & (abs(events["Electron_dz"]) < 0.10 + 0.10 * (abs(events["Electron_eta"]) > 1.479))
            & ((abs(events["Electron_eta"]) < 1.444) | (abs(events["Electron_eta"]) > 1.566))
            & (abs(events["Electron_eta"]) < 2.5)
    )

    # Cut on eta
    eta_split = (
            (np.abs(events["Electron_eta"]) >= eta_min)
            & (np.abs(events["Electron_eta"]) < eta_max)
    )

    # Select based on trigger
    ele = events["Electron_pt"]
    combined_mask = ele_quality_check & eta_split
    evs = ele[combined_mask]
    tr_evs = ele[combined_mask & triggerSingleElectron]

    # Fill histograms
    for pt in ak.flatten(evs):
        hist.Fill(pt)
    for pt in ak.flatten(tr_evs):
        hist.Fill(pt)

    return 0

eta_splits = [[0, 1.0], [1.0, 2.0], [2.0, 3.0]]

# Loop over eta bins and fill histograms
ele_sums(evs, eta_splits[0], eta_ele_sumhist1)
ele_sums(evs, eta_splits[1], eta_ele_sumhist2)
ele_sums(evs, eta_splits[2], eta_ele_sumhist3)

# Set line colors for histograms
eta_ele_sumhist1.SetLineColor(ROOT.kBlack)
eta_ele_sumhist2.SetLineColor(ROOT.kRed)
eta_ele_sumhist3.SetLineColor(ROOT.kBlue)

# Create canvas for sums plot
c1 = ROOT.TCanvas("canvas1", "", 800, 600)

# Create legend
legend = ROOT.TLegend(0.5, 0.6, 0.9, 0.9)
legend.AddEntry(eta_ele_sumhist1, "|#eta|<1.0", "l")
legend.AddEntry(eta_ele_sumhist2, "1.0<|#eta|<2.0", "l")
legend.AddEntry(eta_ele_sumhist3, "2.0<|#eta|<3.0", "l")
legend.AddEntry(ROOT.nullptr, "T = " + temp + "GeV, " + year,"")
legend.AddEntry(ROOT.nullptr, "SUEP decay type: " + decay_type,"")
legend.AddEntry(ROOT.nullptr, "Dark meson mass = " + md,"")
legend.SetTextColor(ROOT.kBlack)
legend.SetTextFont(42)
legend.SetTextSize(0.03)

# Plot sums for eta bins 1, 2, and 3
c1.cd(1)
eta_ele_sumhist1.SetTitle("Number of Electrons in bins of pT;Electron pT [GeV];Number of Electrons")
eta_ele_sumhist1.SetStats(0)
eta_ele_sumhist1.Draw()
eta_ele_sumhist2.Draw("same")
eta_ele_sumhist3.Draw("same")
legend.Draw()

# Save sums plot as PDF
c1.SaveAs(sample_name + "_sums_eta_bins.pdf")

print("Sample " + sample_name + " processed")

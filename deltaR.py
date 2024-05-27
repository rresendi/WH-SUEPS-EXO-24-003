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
    year = "2016 APV conditions"
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
    evs = f['Events'].arrays(['Electron_pt',
                              'Electron_eta',
                              'Electron_phi',
                              'TrigObj_pt',
                              'TrigObj_eta',
                              'TrigObj_phi'])
    return evs

with uproot.open(input_file) as f:
    evs = Events(f)

# Computes deltaR2
def deltaR2(eta1, phi1, eta2, phi2):
    deta = eta1 - eta2
    dphi = phi1 - phi2
    dphi = np.mod(dphi + np.pi, 2*np.pi) - np.pi
    return deta**2 + dphi**2

# Computes deltaR
def deltaR(eta1, phi1, eta2, phi2):
    return np.sqrt(deltaR2(eta1, phi1, eta2, phi2))

# Define Delta R histograms for different eta bins
deltaR_bin_edges = array('d', np.linspace(0, 5, 51))  # 50 bins from 0 to 5

eta1_deltaR_hist = ROOT.TH1D("deltaR_eta_0_0.7", "#Delta R for 0 <= |#eta| < 0.7", len(deltaR_bin_edges)-1, deltaR_bin_edges)
eta2_deltaR_hist = ROOT.TH1D("deltaR_eta_0.7_1.444", "#Delta R for 0.7 <= |#eta| < 1.444", len(deltaR_bin_edges)-1, deltaR_bin_edges)
eta3_deltaR_hist = ROOT.TH1D("deltaR_eta_1.566_2", "#Delta R for 1.566 <= |#eta| < 2", len(deltaR_bin_edges)-1, deltaR_bin_edges)
eta4_deltaR_hist = ROOT.TH1D("deltaR_eta_2_2.5", "#Delta R for 2 <= |#eta| < 2.5", len(deltaR_bin_edges)-1, deltaR_bin_edges)

# Ensure matching shapes and calculate deltaR
electron_trig_pairs = ak.cartesian({"e_eta": evs['Electron_eta'], 
                                    "e_phi": evs['Electron_phi'], 
                                    "t_eta": evs['TrigObj_eta'], 
                                    "t_phi": evs['TrigObj_phi']})

deltaR_values = deltaR(electron_trig_pairs['e_eta'], 
                       electron_trig_pairs['e_phi'], 
                       electron_trig_pairs['t_eta'], 
                       electron_trig_pairs['t_phi'])

# Fill the histograms
for i in range(len(evs['Electron_eta'])):
    for deltaR_value in deltaR_values[i]:
        eta = abs(evs['Electron_eta'][i])
        for eta_value in eta:
            if 0 <= eta_value < 0.7:
                eta1_deltaR_hist.Fill(deltaR_value)
            elif 0.7 <= eta_value < 1.444:
                eta2_deltaR_hist.Fill(deltaR_value)
            elif 1.566 <= eta_value < 2:
                eta3_deltaR_hist.Fill(deltaR_value)
            elif 2 <= eta_value < 2.5:
                eta4_deltaR_hist.Fill(deltaR_value)

# Create a canvas and plot the histograms
c1 = ROOT.TCanvas("canvas", "", 800, 600)
legend = ROOT.TLegend(0.5, 0.1, 0.9, 0.4)

eta1_deltaR_hist.SetLineColor(1)
eta2_deltaR_hist.SetLineColor(2)
eta3_deltaR_hist.SetLineColor(4)
eta4_deltaR_hist.SetLineColor(8)

eta1_deltaR_hist.Draw()
eta2_deltaR_hist.Draw("same")
eta3_deltaR_hist.Draw("same")
eta4_deltaR_hist.Draw("same")

legend.AddEntry(eta1_deltaR_hist, "0 <= |#eta| < 0.7", "l")
legend.AddEntry(eta2_deltaR_hist, "0.7 <= |#eta| < 1.444", "l")
legend.AddEntry(eta3_deltaR_hist, "1.566 <= |#eta| < 2", "l")
legend.AddEntry(eta4_deltaR_hist, "2 <= |#eta| < 2.5", "l")
legend.Draw("same")

c1.Update()
c1.SaveAs(folder + sample_name + "_deltaR_EtaBins.pdf")

# Save histograms to the ROOT file
root_file = ROOT.TFile(output_file, "UPDATE")
root_file.cd()

deltaR_dir = root_file.Get("DeltaR_Histograms")
if not deltaR_dir:
    deltaR_dir = root_file.mkdir("DeltaR_Histograms")
deltaR_dir.cd()

eta1_deltaR_hist.Write()
eta2_deltaR_hist.Write()
eta3_deltaR_hist.Write()
eta4_deltaR_hist.Write()

root_file.Close()

print("sample " + sample_name + " complete")

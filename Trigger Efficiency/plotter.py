import ROOT
import os
import sys
import CMS_lumi
import tdrstyle
from array import array

#####################
### Configuration ###
#####################

# for x label if doing pt
lepton = sys.argv[1]  # Electron or Muon
isdata = sys.argv[2]  # mc or data
eta = sys.argv[3]  # Eta or noEta
era = sys.argv[4]
output = sys.argv[5]  # website?
input_folder = sys.argv[6]  # Input folder containing ROOT files

if lepton == "Muon":
    pt_label = "p_{T}^{#mu} [GeV]"
    phi_label = "#phi^{#mu}"
else:
    pt_label = "p_{T}^{e} [GeV]"
    phi_label = "#phi^{e}"

# Define bin edges for lep1pt as a regular list
lep1pt_bin_edges = array('d', [0, 2, 4, 6, 8, 10, 12,
                             14, 16, 18, 20, 22,
                            24, 26, 28, 30, 32,
                           34, 36, 38, 40, 50,
                          60, 70, 80, 90, 100,
                         120, 140, 160, 180, 200])
                         
histBins = {
    "lep1pt": lep1pt_bin_edges,
    "MET": [30, 0, 300],
    "mT": [15, 0, 150],
    "lep1phi": [35, 0, 3.5],
    "MET phi": [35, 0, 3.5],
    "Boson pT": [80, 0, 800]
}

# TDR
tdrstyle.setTDRStyle()

# Change the CMS_lumi variables (see CMS_lumi.py)
CMS_lumi.cmsText = "CMS"
CMS_lumi.writeExtraText = True

if isdata == "data":
    CMS_lumi.extraText = "Preliminary"
    if era in ["2016", "2016APV"]:
        CMS_lumi.lumi_13TeV = "36.31 fb^{-1}"
    elif era == "2017":
        CMS_lumi.lumi_13TeV = "41.48 fb^{-1}"
    else: #era == 2018
        CMS_lumi.lumi_13TeV = "59.83 fb^{-1}"
else: #mc
    CMS_lumi.extraText = "Simulation Preliminary"

CMS_lumi.lumi_sqrtS = "13 TeV"  # Used with iPeriod = 0, e.g., for simulation-only plots

iPos = 0
if iPos == 0:
    CMS_lumi.relPosX = 0.13

H_ref = 600
W_ref = 800
W = W_ref
H = H_ref

# Define the period for lumi
if isdata == "data":
    iPeriod = 4
else:
    iPeriod = 0  # 13 TeV only

# Function to add histograms from multiple files
def add_histograms(files, hist_name):
    hist = None
    for file_name in files:
        file = ROOT.TFile(file_name, "READ")
        if not file:
            print(f"File {file_name} not found!")
            continue
        h = file.Get(hist_name)
        if not h:
            print(f"Histogram {hist_name} not found in {file_name}!")
            continue
        if hist is None:
            hist = h.Clone()
            hist.SetDirectory(0)  # Disown histogram from TFile
        else:
            hist.Add(h)
        file.Close()
    return hist

# List all ROOT files in the input folder
root_files = [os.path.join(input_folder, f) for f in os.listdir(input_folder) if f.endswith(".root")]

# Iterate through histogram bins and create efficiency plots
for var in histBins:
    num_hist = add_histograms(root_files, var + "_num")
    den_hist = add_histograms(root_files, var + "_den")
    if not num_hist or not den_hist:
        print(f"Histograms for {var} not found in any file!")
        continue

    teff = ROOT.TEfficiency(num_hist, den_hist)
    c = ROOT.TCanvas(var + "C", var + "C", 800, 600)
    c.SetFillColor(0)
    c.SetBorderMode(0)
    c.SetFrameFillStyle(0)
    c.SetFrameBorderMode(0)
    c.SetLeftMargin(0.15)
    c.SetRightMargin(0.04)
    c.SetTopMargin(0.08)
    c.SetBottomMargin(0.18)
    c.SetTickx(1)
    c.SetTicky(1)

    if var == "lep1pt":
        x_label = pt_label
    elif var == "lep1phi":
        x_label = phi_label
    elif var == "Boson pT":
        x_label = "p_{T}^{W} [GeV]"
    else:
        x_label = var

    teff.SetTitle(";%s;%s" % (x_label, "Efficiency"))
    teff.Draw("AP")
    c.Update()
    gr = teff.GetPaintedGraph()
    gr.SetMinimum(0)
    gr.SetMaximum(1.05)
    if var == "lep1pt":
        gr.GetXaxis().SetLimits(0, 200)
    elif var == "Boson pT":
        gr.GetXaxis().SetLimits(0, 800)
    else:
        gr.GetXaxis().SetLimits(histBins[var][1], histBins[var][2])

    # Apply the CMS lumi style
    CMS_lumi.CMS_lumi(c, iPeriod, iPos)

    c.SaveAs(os.path.join(output, f"{lepton}{isdata}{eta}{era}{var}.pdf"))

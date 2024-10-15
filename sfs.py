import ROOT
import os
import sys
import CMS_lumi
import tdrstyle
from array import array

# System Arguments
lepton = sys.argv[1]
isdata = sys.argv[2]
eta = sys.argv[3]
output = sys.argv[4]
data_dir = sys.argv[5]
mc_dir = sys.argv[6]

# Set labels based on the lepton type
if "muon" in lepton.lower():
    pt_label = "Muon p_{T} [GeV]"
    phi_label = "Muon #phi"
else:
    pt_label = "Electron p_{T} [GeV]"
    phi_label = "Electron #phi"

# Define bin edges for 'lep1pt' as a regular list
lep1pt_bin_edges = array('d', [
    0, 2, 4, 6, 8, 10, 12,
    14, 16, 18, 20, 22,
    24, 26, 28, 30, 32,
    34, 36, 38, 40, 50,
    60, 70, 80, 90, 100,
    120, 140, 160, 180, 200
])

# Update histBins dictionary
histBins = {
    "lep1pt": lep1pt_bin_edges,
    "MET": [30, 0, 300],
    "mT": [15, 0, 150],
    "lep1phi": [35, 0, 3.5],
    "MET phi": [35, 0, 3.5]
}

# Set TDR style
tdrstyle.setTDRStyle()

# Change the CMS_lumi variables
CMS_lumi.cmsText = "CMS"
CMS_lumi.writeExtraText = True
CMS_lumi.extraText = "Preliminary"
CMS_lumi.lumi_13TeV = "36.31 fb^{-1}"
CMS_lumi.lumi_sqrtS = "13 TeV"  # Used with iPeriod = 0

iPos = 0
if iPos == 0:
    CMS_lumi.relPosX = 0.13

H_ref = 600
W_ref = 800
W = W_ref
H = H_ref

# Define the period for lumi
iPeriod = 0  # 13 TeV only

# loading and combining histograms from all ROOT files in a directory
def get_combined_histogram(directory, hist_name):
    combined_hist = None
    for filename in os.listdir(directory):
        if filename.endswith(".root"):
            filepath = os.path.join(directory, filename)
            file = ROOT.TFile(filepath, "READ")
            hist = file.Get(hist_name)

            # Check if the histogram is found and has entries
            if hist:
                print(f"Found histogram '{hist_name}' in {filename}")
                if hist.GetEntries() == 0:
                    print(f"Warning: Histogram '{hist_name}' in {filename} has no entries.")
                if combined_hist is None:
                    combined_hist = hist.Clone()
                    combined_hist.SetDirectory(0)
                else:
                    combined_hist.Add(hist)
            else:
                print(f"Warning: Could not find histogram '{hist_name}' in {filename}")
            file.Close()

    # If no histograms were found
    if combined_hist is None:
        print(f"Error: No histograms named '{hist_name}' found in directory {directory}")
    return combined_hist

# Create output directory if it doesn't exist
if not os.path.exists(output):
    os.makedirs(output)

# List of variables to plot
variables = histBins.keys()

# Iterate through histogram bins and create efficiency plots
for var in variables:
    # Combine numerator and denominator histograms from all data and MC files
    data_num_hist = get_combined_histogram(data_dir, var + "_num")
    data_den_hist = get_combined_histogram(data_dir, var + "_den")
    mc_num_hist = get_combined_histogram(mc_dir, var + "_num")
    mc_den_hist = get_combined_histogram(mc_dir, var + "_den")

    # Check if the histograms exist and have entries
    if not data_num_hist or not data_den_hist or not mc_num_hist or not mc_den_hist:
        print(f"Error: Could not find necessary histograms for variable '{var}' in data or MC directories, or histograms are empty.")
        continue

    # Create TEfficiency objects for data and MC
    data_eff = ROOT.TEfficiency(data_num_hist, data_den_hist)
    mc_eff = ROOT.TEfficiency(mc_num_hist, mc_den_hist)

    # Create canvas
    c = ROOT.TCanvas(var + "_c", var + "_c", W, H)
    c.SetFillColor(0)
    c.SetBorderMode(0)
    c.SetFrameFillStyle(0)
    c.SetFrameBorderMode(0)
    c.SetLeftMargin(0.15)
    c.SetRightMargin(0.04)
    c.SetTopMargin(0.08)
    c.SetBottomMargin(0.18)
    c.SetTickx(0)
    c.SetTicky(0)

    # Set x-axis label
    if var == "lep1pt":
        x_label = pt_label
    elif var == "lep1phi":
        x_label = phi_label
    else:
        x_label = var

    # Draw data efficiency
    data_eff.SetTitle(f";{x_label};Efficiency")
    data_eff.SetMarkerStyle(20)
    data_eff.SetMarkerColor(ROOT.kBlack)
    data_eff.SetLineColor(ROOT.kBlack)
    data_eff.Draw("AP")

    # Draw MC efficiency on the same canvas
    mc_eff.SetMarkerStyle(21)
    mc_eff.SetMarkerColor(ROOT.kRed)
    mc_eff.SetLineColor(ROOT.kRed)
    mc_eff.Draw("P same")

    # Set axis ranges
    c.Update()
    gr = data_eff.GetPaintedGraph()
    gr.SetMinimum(0)
    gr.SetMaximum(1.05)
    if var == "lep1pt":
        gr.GetXaxis().SetLimits(0, 200)
    else:
        gr.GetXaxis().SetLimits(histBins[var][1], histBins[var][2])

    # Add legend
    legend = ROOT.TLegend(0.65, 0.15, 0.85, 0.30)
    legend.AddEntry(data_eff, "Data", "lep")
    legend.AddEntry(mc_eff, "MC", "lep")
    legend.Draw()

    # Apply the CMS lumi style
    CMS_lumi.CMS_lumi(c, iPeriod, iPos)

    # Save the canvas
    c.SaveAs(os.path.join(output, f"{var}_efficiency.pdf"))

    # Now create the ratio of efficiencies
    # Create histograms from TEfficiency objects
    data_eff_hist = data_eff.GetCopyTotalHisto()
    data_eff_hist.Reset()
    mc_eff_hist = mc_eff.GetCopyTotalHisto()
    mc_eff_hist.Reset()

    # Fill histograms with efficiencies and errors
    for i in range(1, data_eff_hist.GetNbinsX() + 1):
        if data_eff.GetTotalHistogram().GetBinContent(i) > 0:
            eff_data = data_eff.GetEfficiency(i)
            eff_data_err_low = data_eff.GetEfficiencyErrorLow(i)
            eff_data_err_up = data_eff.GetEfficiencyErrorUp(i)
            data_eff_hist.SetBinContent(i, eff_data)
            data_eff_hist.SetBinError(i, max(eff_data_err_low, eff_data_err_up))
        else:
            data_eff_hist.SetBinContent(i, 0)
            data_eff_hist.SetBinError(i, 0)

        if mc_eff.GetTotalHistogram().GetBinContent(i) > 0:
            eff_mc = mc_eff.GetEfficiency(i)
            eff_mc_err_low = mc_eff.GetEfficiencyErrorLow(i)
            eff_mc_err_up = mc_eff.GetEfficiencyErrorUp(i)
            mc_eff_hist.SetBinContent(i, eff_mc)
            mc_eff_hist.SetBinError(i, max(eff_mc_err_low, eff_mc_err_up))
        else:
            mc_eff_hist.SetBinContent(i, 0)
            mc_eff_hist.SetBinError(i, 0)

    # Calculate ratio histogram
    ratio_hist = data_eff_hist.Clone(var + "_efficiency_ratio")
    ratio_hist.Divide(mc_eff_hist)

    # Create canvas for ratio
    c_ratio = ROOT.TCanvas(var + "_ratio_c", var + "_ratio_c", W, H)
    c_ratio.SetFillColor(0)
    c_ratio.SetBorderMode(0)
    c_ratio.SetFrameFillStyle(0)
    c_ratio.SetFrameBorderMode(0)
    c_ratio.SetLeftMargin(0.15)
    c_ratio.SetRightMargin(0.04)
    c_ratio.SetTopMargin(0.08)
    c_ratio.SetBottomMargin(0.18)
    c_ratio.SetTickx(0)
    c_ratio.SetTicky(0)

    # Set axis titles and ranges
    ratio_hist.SetTitle(f";{x_label};Data/MC Efficiency Ratio")
    ratio_hist.SetMinimum(0.0)
    ratio_hist.SetMaximum(2.0)
    ratio_hist.SetMarkerStyle(20)
    ratio_hist.SetMarkerColor(ROOT.kBlue)
    ratio_hist.SetLineColor(ROOT.kBlue)

    # Draw ratio histogram
    ratio_hist.Draw("EP")

    # Apply CMS lumi style
    CMS_lumi.CMS_lumi(c_ratio, iPeriod, iPos)

    # Save the ratio plot
    c_ratio.SaveAs(os.path.join(output, f"{var}_efficiency_ratio.pdf"))

print("Processing complete. Plots saved in:", output)

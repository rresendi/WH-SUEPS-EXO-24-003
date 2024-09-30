import ROOT
import os
import sys
import CMS_lumi
import tdrstyle
from array import array

#####################
### Configuration ###
#####################

# Get input parameters from the command line
lepton = sys.argv[1]  # Electron or Muon
isdata = sys.argv[2]  # mc or data
output = sys.argv[3]  # Output directory for plot
data_dir = sys.argv[4]  # Data directory containing ROOT files
mc_dir = sys.argv[5]  # MC directory containing ROOT files

# Define bin edges for lep1pt as provided
lep1pt_bin_edges = array('d', [0, 2, 4, 6, 8, 10, 12,
                               14, 16, 18, 20, 22,
                               24, 26, 28, 30, 32,
                               34, 36, 38, 40, 50,
                               60, 70, 80, 90, 100,
                               120, 140, 160, 180, 200])

# Set CMS TDR Style
tdrstyle.setTDRStyle()

# CMS lumi configuration
CMS_lumi.cmsText = "CMS"
CMS_lumi.writeExtraText = True
CMS_lumi.extraText = "Preliminary"
CMS_lumi.lumi_13TeV = "36.31 fb^{-1}"
CMS_lumi.lumi_sqrtS = "13 TeV"  # For simulation-only plots

iPos = 0
if iPos == 0:
    CMS_lumi.relPosX = 0.13

# Canvas dimensions
H_ref = 600
W_ref = 800
W = W_ref
H = H_ref

# Define the period for lumi (13 TeV)
iPeriod = 0  # Change this to use other datasets, e.g., 7 or 8 TeV

##########################
### Helper Functions ###
##########################

# Function to load and combine histograms from all ROOT files in a directory
def get_combined_histogram(directory, hist_name):
    combined_hist = None
    for filename in os.listdir(directory):
        if filename.endswith(".root"):
            filepath = os.path.join(directory, filename)
            file = ROOT.TFile(filepath, "READ")
            hist = file.Get(hist_name)

            # Debugging: Check if the histogram is found and has entries
            if hist:
                print(f"Found histogram '{hist_name}' in {filename}")
                if hist.GetEntries() == 0:
                    print(f"Warning: Histogram '{hist_name}' in {filename} has no entries.")
                if combined_hist is None:
                    combined_hist = hist.Clone()
                    combined_hist.SetDirectory(0)  # Detach from file to keep it in memory
                else:
                    combined_hist.Add(hist)  # Combine histograms
            else:
                print(f"Warning: Could not find histogram '{hist_name}' in {filename}")
            file.Close()

    # If no histograms were found, print an error
    if combined_hist is None:
        print(f"Error: No histograms found in directory {directory}")
    return combined_hist

##########################
### Plotting Procedure ###
##########################

# Combine histograms from all data and MC files
data_hist = get_combined_histogram(data_dir, "lep1pt_num")
mc_hist = get_combined_histogram(mc_dir, "lep1pt_num")

# Check if the histograms exist and have entries
if not data_hist or not mc_hist:
    print("Error: Could not find lep1pt histograms in either data or MC directory, or histograms are empty.")
    sys.exit(1)

# Debugging: Check number of entries in histograms
print(f"Data histogram 'lep1pt_num' has {data_hist.GetEntries()} entries.")
print(f"MC histogram 'lep1pt_num' has {mc_hist.GetEntries()} entries.")

# Create the ratio histogram (h1/h2 with uncertainties automatically propagated)
ratio_hist = data_hist.Clone("lep1pt_ratio")
ratio_hist.Divide(mc_hist)  # This handles both the division and the uncertainty propagation

# Debugging: Print bin contents for the ratio histogram
for bin_idx in range(1, ratio_hist.GetNbinsX() + 1):
    print(f"Bin {bin_idx}: Data/MC Ratio = {ratio_hist.GetBinContent(bin_idx)}")

# Create a canvas for the plot
c = ROOT.TCanvas("lep1ptC", "lep1ptC", W, H)
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

# Apply binning to the ratio histogram
ratio_hist.SetTitle(";Lepton pT [GeV];Data/MC Ratio")
ratio_hist.SetMinimum(0)
ratio_hist.SetMaximum(15)  # Adjust as needed to show the correct range

# Draw the ratio plot with error bars ('E' option for errors)
ratio_hist.Draw("E")  # 'E' includes error bars for uncertainties

# Apply the CMS lumi style
CMS_lumi.CMS_lumi(c, iPeriod, iPos)

# Save the plot as a PDF
output_filename = f"{output}/lep1pt_ratio.pdf"
c.SaveAs(output_filename)
print(f"Saved plot for lep1pt ratio to {output_filename}")

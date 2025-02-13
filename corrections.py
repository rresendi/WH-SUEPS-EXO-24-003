import ROOT
import os
import sys
import json
from array import array

# System Arguments
lepton = sys.argv[1]
era = sys.argv[2]
output = sys.argv[3]
data_dir = sys.argv[4]
mc_dir = sys.argv[5]

# Set labels based on the lepton type
if lepton == "Muon":
    pt_label = "p_{T}^{#mu} [GeV]"
    phi_label = "#phi^{#mu}"
else:
    pt_label = "p_{T}^{e} [GeV]"
    phi_label = "#phi^{e}"

# Define bin edges for 'lep1pt'
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

# Loading and combining histograms from all ROOT files in a directory
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

    if combined_hist is None:
        print(f"Error: No histograms named '{hist_name}' found in directory {directory}")
    return combined_hist

# Create output directory if it doesn't exist
if not os.path.exists(output):
    os.makedirs(output)

# List of variables to process
variables = histBins.keys()

# Dictionary to store scale factors for all variables
scale_factors = {}

# Iterate through histogram bins and compute scale factors
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

    # Compute scale factors and store bin edges
    scale_factors[var] = []
    n_bins = data_eff.GetTotalHistogram().GetNbinsX()

    for i in range(1, n_bins + 1):
        # Bin edges
        low_edge = data_eff.GetTotalHistogram().GetXaxis().GetBinLowEdge(i)
        high_edge = data_eff.GetTotalHistogram().GetXaxis().GetBinUpEdge(i)

        # Efficiencies
        eff_data = data_eff.GetEfficiency(i)
        eff_mc = mc_eff.GetEfficiency(i)

        # Compute scale factor (handle division by zero)
        sf = eff_data / eff_mc if eff_mc > 0 else 0.0

        # Store bin edges and scale factor
        scale_factors[var].append({
            "bin_range": [low_edge, high_edge],
            "scale_factor": sf
        })

# Write scale factors to a JSON file
output_json_file = os.path.join(output, f"{lepton}{era}_scale_factors.json")
with open(output_json_file, 'w') as f_out:
    json.dump(scale_factors, f_out, indent=2)

print(f"Scale factors written to {output_json_file}")

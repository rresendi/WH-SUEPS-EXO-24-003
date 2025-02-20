import ROOT
import os
import sys
import json
from array import array

# System Arguments
lepton = sys.argv[1]
era = sys.argv[2]
output = sys.argv[3]
data_ref_dir = sys.argv[4]
mc_dir = sys.argv[5]
data_noref_dir = sys.argv[6]

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

# Function to load and combine histograms from all ROOT files in a directory
def get_combined_histogram(directory, hist_name):
    combined_hist = None
    for filename in os.listdir(directory):
        if filename.endswith(".root"):
            filepath = os.path.join(directory, filename)
            file = ROOT.TFile(filepath, "READ")
            hist = file.Get(hist_name)

            if hist:
                if combined_hist is None:
                    combined_hist = hist.Clone()
                    combined_hist.SetDirectory(0)
                else:
                    combined_hist.Add(hist)
            file.Close()

    return combined_hist

# Create output directory if it doesn't exist
if not os.path.exists(output):
    os.makedirs(output)

# Dictionary to store scale factors for all variables
scale_factors = {}

# Iterate through histogram bins and compute scale factors
for var in histBins.keys():
    # Load histograms for Data (Reference), MC, and Data (No-Reference)
    data_ref_num = get_combined_histogram(data_ref_dir, var + "_num")
    data_ref_den = get_combined_histogram(data_ref_dir, var + "_den")
    mc_num = get_combined_histogram(mc_dir, var + "_num")
    mc_den = get_combined_histogram(mc_dir, var + "_den")
    data_noref_num = get_combined_histogram(data_noref_dir, var + "_num")
    data_noref_den = get_combined_histogram(data_noref_dir, var + "_den")

    # Check if all required histograms exist
    if not data_ref_num or not data_ref_den or not mc_num or not mc_den or not data_noref_num or not data_noref_den:
        print(f"Error: Missing histograms for variable '{var}'. Skipping...")
        continue

    # Create TEfficiency objects
    eff_data_ref = ROOT.TEfficiency(data_ref_num, data_ref_den)
    eff_mc = ROOT.TEfficiency(mc_num, mc_den)
    eff_data_noref = ROOT.TEfficiency(data_noref_num, data_noref_den)

    # Compute scale factors
    scale_factors[var] = []
    n_bins = eff_data_ref.GetTotalHistogram().GetNbinsX()

    for i in range(1, n_bins + 1):
        # Bin edges
        low_edge = eff_data_ref.GetTotalHistogram().GetXaxis().GetBinLowEdge(i)
        high_edge = eff_data_ref.GetTotalHistogram().GetXaxis().GetBinUpEdge(i)

        # Efficiencies
        eff_ref = eff_data_ref.GetEfficiency(i)
        eff_mc_val = eff_mc.GetEfficiency(i)
        eff_noref = eff_data_noref.GetEfficiency(i)

        # Compute Clopper-Pearson statistical uncertainties
        err_ref_low = eff_data_ref.GetEfficiencyErrorLow(i)
        err_ref_up = eff_data_ref.GetEfficiencyErrorUp(i)
        err_mc_low = eff_mc.GetEfficiencyErrorLow(i)
        err_mc_up = eff_mc.GetEfficiencyErrorUp(i)
        err_noref_low = eff_data_noref.GetEfficiencyErrorLow(i)
        err_noref_up = eff_data_noref.GetEfficiencyErrorUp(i)

        # Compute scale factor
        sf = eff_ref / eff_mc_val if eff_mc_val > 0 else 0.0

        # Compute reference/no-reference efficiency ratio
        ref_noref_ratio = eff_ref / eff_noref if eff_noref > 0 else 1.0

        # Error propagation for reference/no-reference ratio
        err_ratio_low = ref_noref_ratio * ((err_ref_low / eff_ref) ** 2 + (err_noref_up / eff_noref) ** 2) ** 0.5 if eff_ref > 0 and eff_noref > 0 else 0.0
        err_ratio_up = ref_noref_ratio * ((err_ref_up / eff_ref) ** 2 + (err_noref_low / eff_noref) ** 2) ** 0.5 if eff_ref > 0 and eff_noref > 0 else 0.0

        # Error propagation for scale factor (MC and Data Ref errors)
        sf_err_low = sf * ((err_ref_low / eff_ref) ** 2 + (err_mc_up / eff_mc_val) ** 2) ** 0.5 if eff_ref > 0 and eff_mc_val > 0 else 0.0
        sf_err_up = sf * ((err_ref_up / eff_ref) ** 2 + (err_mc_low / eff_mc_val) ** 2) ** 0.5 if eff_ref > 0 and eff_mc_val > 0 else 0.0

        # Quadratically combine the statistical and systematic uncertainties
        total_err_low = (sf_err_low**2 + err_ratio_low**2) ** 0.5
        total_err_up = (sf_err_up**2 + err_ratio_up**2) ** 0.5

        # Store values in JSON output
        scale_factors[var].append({
            "bin_range": [low_edge, high_edge],
            # "scale_factor": sf,
            # "total_uncertainty_low": total_err_low,
            # "total_uncertainty_up": total_err_up,
            "systematic uncertainty_low": err_ratio_low,
            "systematic uncertainty_up": err_ratio_up
        })

# Write scale factors to a JSON file
output_json_file = os.path.join(output, f"{lepton}{era}_sys_uncertainties.json")
with open(output_json_file, 'w') as f_out:
    json.dump(scale_factors, f_out, indent=2)

print(f"Scale factors written to {output_json_file}")

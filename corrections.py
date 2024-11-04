import ROOT
import numpy as np
import sys
import json
import correctionlib.schemav2 as cs

# Define the lep1pt bins
lep1pt_bins = [0, 2, 4, 6, 8, 10, 12,
               14, 16, 18, 20, 22,
               24, 26, 28, 30, 32,
               34, 36, 38, 40, 50,
               60, 70, 80, 90, 100,
               120, 140, 160, 180, 200]

# Check for input ROOT file argument
if len(sys.argv) != 2:
    print("Usage: python script.py <input_root_file>")
    sys.exit(1)

input_root_file = sys.argv[1]

# Open the ROOT file
f = ROOT.TFile.Open(input_root_file)
if not f or f.IsZombie():
    print("Error opening file:", input_root_file)
    sys.exit(1)

# Get the numerator and denominator histograms
num_hist = f.Get('lep1pt_num')
den_hist = f.Get('lep1pt_den')

# Check if histograms are found
if not num_hist:
    print("Error: 'lep1pt_num' histogram not found in the ROOT file.")
    sys.exit(1)
if not den_hist:
    print("Error: 'lep1pt_den' histogram not found in the ROOT file.")
    sys.exit(1)

# Compute scale factors per bin
scale_factors = []
num_bins = num_hist.GetNbinsX()
for i in range(1, num_bins + 1):
    num_content = num_hist.GetBinContent(i)
    den_content = den_hist.GetBinContent(i)
    if den_content != 0:
        sf = num_content / den_content
    else:
        sf = 0.0  # Handle zero denominator
    scale_factors.append(sf)

# Ensure the number of scale factors matches the number of bins
if len(scale_factors) != len(lep1pt_bins) - 1:
    raise ValueError("Number of scale factors does not match number of bins")

# Create the binning evaluator
binned_evaluator = cs.Binning(
    nodetype='Binning',
    input='lep1pt',
    edges=lep1pt_bins,
    content=scale_factors,
    flow='clamp'
)

# Create the correction object
correction = cs.Correction(
    version=1,
    name='lep1pt_scale_factor',
    description='Scale factors as a function of lep1pt',
    outputs=[cs.Output(name='weight', type='real')],
    inputs=[cs.Input(name='lep1pt', type='real')],
    data=binned_evaluator
)

# Create the correction set
correction_set = cs.CorrectionSet(
    schema_version=2,
    corrections=[correction]
)

# Write the correction set to a JSON file
output_json_file = input_root_file.replace('.root', '_corrections.json')
with open(output_json_file, 'w') as f_out:
    json.dump(correction_set.dict(), f_out, indent=2)

print(f"Corrections written to {output_json_file}")

import argparse
import re
import matplotlib.pyplot as plt
import numpy as np

def readSAFFile(safContent, section):
    """
    Parse the SAF file content for a given histogram section.
    """
    pat = rf'<Histo>[\s\S]*?<Description>\s*"{section}"[\s\S]*?<Data>([\s\S]*?)</Data>'
    m = re.search(pat, safContent)
    if not m:
        return []
    lines = m.group(1).strip().split('\n')
    if len(lines) < 3:
        return []
    values = [float(line.split()[0]) for line in lines]
    # Adjust first and last bins as in the original procedure.
    values[1] += values[0]
    values[-2] += values[-1]
    return values[1:-1]

def plotKinematic(bin_centers, yields, uncertainties, kinematic):
    # Check that the lengths of bin_centers, yields, and uncertainties are the same
    if len(bin_centers) != len(yields) or len(bin_centers) != len(uncertainties):
        print(f"Error: Mismatch in lengths for {kinematic}.")
        return

    # Plot the kinematic variable (raw counts with uncertainties)
    plt.errorbar(bin_centers, yields, yerr=uncertainties, fmt='o', label=kinematic, color='navy', ecolor='navy')
    plt.xlabel(f'{kinematic}')
    plt.ylabel('Number of Events')
    plt.title(f'{kinematic} Plot')
    plt.legend()
    plt.tight_layout()

    # Dynamically adjust the x-axis range based on the kinematic value
    plt.xlim(min(bin_centers) - 0.5, max(bin_centers) + 0.5)
    outname = f"/eos/user/g/gdecastr/www/WH_Reinterpretation/FirstPlots/Hadronic_2_2_OGCard/{kinematic}.png"
    plt.savefig(outname)
    plt.close()
    print(f"Saved plot: {outname}")

def main():
    parser = argparse.ArgumentParser(description="Plot kinematics (raw number of events with uncertainties)")
    parser.add_argument("--saf", required=True, help="Path to the SAF file")
    args = parser.parse_args()
    
    saf_file = args.saf

    with open(saf_file, "r") as f:
        saf_content = f.read()

    # List of kinematics (SAF section names)
    kinematics = [
        "ABCD_A",
        "ABCD_B",
        "ABCD_C",
        "ABCD_D",
        "ABCD_E",
        "ABCD_F0",
        "ABCD_F1",
        "ABCD_F2",
        "ABCD_F3",
        "ABCD_F4",
        "ABCD_H",
        "ABCD_SR0",
        "ABCD_SR1",
        "ABCD_SR2",
        "ABCD_SR3",
        "ABCD_SR4",
        "wMass",
        "wPt",
        "wEta",
        "wPhi",
        "lepPt",
        "lepEta",
        "lepPhi",
        "looseLep",
        "tightLep",
        "muPt",
        "muEta",
        "muPhi",
        "looseMu",
        "tightMu",
        "elePt",
        "eleEta",
        "elePhi",
        "looseEle",
        "tightEle",
        "ak4Pt",
        "ak4Eta", 
        "ak4Phi",
        "ak4NTracks",
        "ak15Pt",
        "ak15Eta",
        "ak15Phi",
        "ak15NTracks",
        "ak15Mass",
        "NJets",
        "metPt",
        "metEta",
        "metPhi",
        "labSphericity",
        "boostedSphericity"
    ]
    
    # Binning for each kinematic
    binning_dict = {
        "ABCD_A": (10, 10.0, 20.0),
        "ABCD_B": (10, 20.0, 30.0),
        "ABCD_C": (170, 30.0, 200.0),
        "ABCD_D": (10, 10.0, 20.0),
        "ABCD_E": (10, 20.0, 30.0),
        "ABCD_F0": (10, 30.0, 40.0),
        "ABCD_F1": (10, 40.0, 50.0),
        "ABCD_F2": (10, 50.0, 60.0),
        "ABCD_F3": (20, 60.0, 80.0),
        "ABCD_F4": (120, 80.0, 200.0),
        "ABCD_H": (10, 20.0, 30.0),
        "ABCD_SR0": (10, 30.0, 40.0),
        "ABCD_SR1": (10, 40.0, 50.0),
        "ABCD_SR2": (10, 50.0, 60.0),
        "ABCD_SR3": (20, 60.0, 80.0),
        "ABCD_SR4": (20, 80.0, 100.0),
        "wMass": (60, 0, 140),
        "wPt": (300, 0, 300),
        "wEta": (40, -3.14, 3.14),
        "wPhi": (40, -3.14, 3.14),
        "lepPt": (300, 0, 300),
        "lepEta": (40, -3.14, 3.14),
        "lepPhi": (40, -3.14, 3.14),
        "looseLep": (10, 0, 10),
        "tightLep": (10, 0, 10),
        "muPt": (300, 0, 300),
        "muEta": (40, -3.14, 3.14),
        "muPhi": (40, -3.14, 3.14),
        "looseMu": (10, 0, 10),
        "tightMu": (10, 0, 10),
        "elePt": (300, 0, 300),
        "eleEta": (40, -3.14, 3.14),
        "elePhi": (40, -3.14, 3.14),
        "looseEle": (10, 0, 10),
        "tightEle": (10, 0, 10),
        "ak4Pt": (300, 0, 300),
        "ak4Eta": (40, -3.14, 3.14),
        "ak4Phi": (40, -3.14, 3.14),
        "ak4NTracks": (100, 0, 100),
        "ak15Pt": (300, 0, 300),
        "ak15Eta": (40, -3.14, 3.14),
        "ak15Phi": (40, -3.14, 3.14),
        "ak15NTracks": (100, 0, 100),
        "ak15Mass": (400, 0, 400),
        "metPt": (300, 0, 300),
        "metEta": (40, -3.14, 3.14),
        "metPhi": (40, -3.14, 3.14),
        "labSphericity": (50, 0, 1),
        "boostedSphericity": (50, 0, 1),
        "NJets": (10, 0, 10),  
    }

    # Loop over each kinematic section and create a plot.
    for kin in kinematics:
        data = readSAFFile(saf_content, kin)
        if not data:
            print(f"Warning: No data found for {kin}")
            continue
        
        # Calculate statistical uncertainties (assuming Poisson errors).
        uncertainties = np.sqrt(np.array(data))
        
        # Use raw counts as yields (no normalization)
        yields = np.array(data)

        # Fetch the appropriate binning for the current kinematic
        if kin in binning_dict:
            num_bins, bin_min, bin_max = binning_dict[kin]
            
            # Calculate bin edges and centers
            bin_edges = np.linspace(bin_min, bin_max, num_bins + 1)
            bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2

            # Plot the data with error bars
            plotKinematic(bin_centers, yields, uncertainties, kin)
        else:
            print(f"Warning: No binning defined for {kin}, skipping...")

if __name__ == "__main__":
    main()

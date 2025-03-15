#!/usr/bin/env python3
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

def plotKinematic(xBins, yields, uncertainties, kinematic):
    plt.errorbar(xBins, yields, yerr=uncertainties, fmt='o', label=kinematic)
    plt.xlabel('Bin')
    plt.ylabel('Normalized Yield')
    plt.title(kinematic)
    plt.legend()
    plt.tight_layout()
    outname = f"{kinematic}.pdf"
    plt.savefig(outname)
    plt.close()
    print(f"Saved plot: {outname}")

def main():
    parser = argparse.ArgumentParser(description="Plot normalized kinematics")
    parser.add_argument("--saf", required=True, help="Path to the SAF file")
    parser.add_argument("--lumi", type=float, required=True, help="Luminosity")
    parser.add_argument("--xsec", type=float, required=True, help="Cross section")
    parser.add_argument("--total_gen", type=float, required=True, help="Total generated events")
    args = parser.parse_args()
    
    saf_file = args.saf
    lumi = args.lumi
    xsec = args.xsec
    total_gen = args.total_gen
    
    norm = lumi * xsec

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
        "wmass",
        "wpt",
        "weta",
        "wphi",
        "lep1pt",
        "lep1eta",
        "lep1phi",
        "ak41pt",
        "ak41eta",
        "ak41phi",
        "ak41ntracks",
        "ak151pt",
        "ak151eta",
        "ak151phi",
        "ak151ntracks",
        "ak151mass",
        "muo1pt",
        "muo1eta",
        "muo1phi",
        "ele1pt",
        "ele1eta",
        "ele1phi",
        "NJets",
        "metpt",
        "meteta",
        "metphi",
        "sphericity"
    ]

    with open(saf_file, "r") as f:
        saf_content = f.read()

    # Loop over each kinematic section and create a plot.
    for kin in kinematics:
        data = readSAFFile(saf_content, kin)
        if not data:
            print(f"Warning: No data found for {kin}")
            continue
        
        # Calculate statistical uncertainties (assuming Poisson errors).
        uncertainties = np.sqrt(np.array(data))
        normed_yields = np.array(data) * norm / total_gen
        normed_uncertainties = uncertainties * norm / total_gen

        xBins = np.arange(1, len(normed_yields) + 1)
        plotKinematic(xBins, normed_yields, normed_uncertainties, kin)

if __name__ == "__main__":
    main()

import uproot
import matplotlib.pyplot as plt
import numpy as np
import re
import argparse
import os
import gc

# Pair each (ROOT, SAF) histo with binning (nbins, xmin, xmax)
paired_histos = [
    ("electron_eta_SR", "eleEta", (100, -5.0, 5.0)),
    ("electron_phi_SR", "elePhi", (60, -3.2, 3.2)),
    ("electron_pt_SR", "elePt", (500, 0, 500)),
    ("SUEP_pt_SR", "ak15Pt", (100, 0, 1000)),
    ("lepton_eta_SR", "lepEta", (100, -5.0, 5.0)),
    ("lepton_phi_SR", "lepPhi", (60, -3.2, 3.2)),
    ("lepton_pt_SR", "lepPt", (1000, 0, 1000)),
    ("muon_eta_SR", "muEta", (100, -5.0, 5.0)),
    ("muon_phi_SR", "muPhi", (60, -3.2, 3.2)),
    ("muon_pt_SR", "muPt", (500, 0, 500)),
    ("ngood_ak4jets_SR", "NJets", (20, 0, 20)),
    ("SUEP_eta_SR", "ak15Eta", (100, -5.0, 5.0)),
    ("SUEP_mass_SR", "ak15Mass", (150, 0, 2000)),
    ("SUEP_nconst_SR", "ak15NTracks", (200, 0, 200)),
    ("SUEP_phi_SR", "ak15Phi", (100, -6.5, 6.5)),
    ("SUEP_S1_SR", "boostedSphericity", (100, 0, 1)),
    ("W_phi_SR", "wPhi", (60, -3.2, 3.2)),
    ("W_pt_SR", "wPt", (200, 0, 2000)),
]

def readSAFFile(safContent, section):
    pat = rf'<Histo>[\s\S]*?<Description>\s*"{section}"[\s\S]*?<Data>([\s\S]*?)</Data>'
    m = re.search(pat, safContent)
    if not m:
        return np.array([])
    lines = m.group(1).strip().split('\n')
    if len(lines) < 3:
        return np.array([])
    values = [float(line.split()[0]) for line in lines]
    values[1] += values[0]
    values[-2] += values[-1]
    return np.array(values[1:-1])

def get_hist_from_root(root_file, hist_name):
    for k in root_file.keys():
        if hist_name == k.split(";")[0]:
            return root_file[k].to_hist()
    return None

def main():
    parser = argparse.ArgumentParser(description="Overlay SAF vs ROOT histograms with ratio")
    parser.add_argument("--saf", required=True, help="SAF file")
    parser.add_argument("--root", required=True, help="ROOT file")
    parser.add_argument("--outdir", default="ratio_plots", help="Output directory")
    parser.add_argument("--lumi", type=float, required=True, help="Luminosity")
    parser.add_argument("--xsec", type=float, required=True, help="Cross section")
    parser.add_argument("--total_gen", type=float, required=True, help="Total generated events")
    args = parser.parse_args()

    norm_factor = args.lumi * args.xsec / args.total_gen

    os.makedirs(args.outdir, exist_ok=True)

    with open(args.saf, "r") as f:
        saf_content = f.read()

    root_file = uproot.open(args.root)

    for root_name, saf_name, (nbins, xmin, xmax) in paired_histos:
        saf_data = readSAFFile(saf_content, saf_name)
        if saf_data.size == 0:
            print(f"Missing SAF data for {saf_name}")
            continue

        root_hist = get_hist_from_root(root_file, root_name)
        if root_hist is None:
            print(f"Missing ROOT histogram {root_name}")
            continue

        root_vals = root_hist.values()
        if len(root_vals) != nbins:
            print(f"Bin mismatch in {root_name} vs {saf_name}, skipping.")
            continue

        bin_edges = np.linspace(xmin, xmax, nbins + 1)
        bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2

        # Normalize SAF data with scaling
        saf_yields = saf_data * norm_factor
        saf_uncertainties = np.sqrt(saf_data) * norm_factor

        bin_widths = np.diff(bin_edges)
        total_yield = np.sum(saf_yields * bin_widths)
        if total_yield > 0:
            saf_yields /= total_yield
            saf_uncertainties /= total_yield
    
        fig, (ax1, ax2) = plt.subplots(2, 1, gridspec_kw={'height_ratios': [3, 1]}, sharex=True, figsize=(7, 6))

        ax1.hist(bin_centers, bins=bin_edges, weights=root_vals, histtype='step', color='darkgreen', label=f'CMSSW: {root_name}', density = True)
        ax1.errorbar(bin_centers, saf_yields, yerr=saf_uncertainties, fmt='o', color='navy', label=f"Delphes: {saf_name}")
        ax1.set_ylabel('Events')
        ax1.legend()
        ax1.set_title(f"{saf_name} vs {root_name}")

        # Ratio (SAF / ROOT)
        with np.errstate(divide='ignore', invalid='ignore'):
            ratio = np.true_divide(saf_yields, root_vals)
            ratio[~np.isfinite(ratio)] = 0

        ax2.step(bin_centers, ratio, where='mid', color='black')
        ax2.axhline(1, color='gray', linestyle='--')
        ax2.set_ylabel('Ratio')
        ax2.set_xlabel(saf_name)

        plt.tight_layout()
        out_path = os.path.join(args.outdir, f"{saf_name}_vs_{root_name}.pdf")
        plt.savefig(out_path)
        plt.close()
        print(f"Saved plot: {out_path}")

if __name__ == "__main__":
    main()

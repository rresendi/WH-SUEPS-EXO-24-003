import uproot
import matplotlib.pyplot as plt
import numpy as np
import re
import argparse
import os
import gc
import ROOT as rpy

# Pair each (ROOT, SAF) histo with binning (nbins, xmin, xmax, rebin)
paired_histos = [
    ("lepton_pt_SR", "lepPt", (370, 30, 400, 5)),
    ("lepton_eta_SR", "lepEta", (50, -2.5, 2.5, 3)),
    ("lepton_phi_SR", "lepPhi", (60, -3.2, 3.2, 3)),

    ("electron_pt_SR", "elePt", (370, 30, 400, 5)),
    ("electron_eta_SR", "eleEta", (50, -2.5, 2.5, 3)),
    ("electron_phi_SR", "elePhi", (60, -3.2, 3.2, 3)),

    ("muon_pt_SR", "muPt", (370, 30, 400, 5)),
    ("muon_eta_SR", "muEta", (50, -2.5, 2.5, 3)),
    ("muon_phi_SR", "muPhi", (60, -3.2, 3.2, 3)),

    ("ngood_ak4jets_SR", "NJets", (20, 0, 20, 1)),

    ("SUEP_pt_SR", "ak15Pt", (44, 60, 500, 5)),
    ("SUEP_eta_SR", "ak15Eta", (50, -2.5, 2.5, 3)),
    ("SUEP_phi_SR", "ak15Phi", (50, -3.25, 3.25, 3)),
    ("SUEP_mass_SR", "ak15Mass", (30, 0, 400, 5)),
    ("SUEP_nconst_SR", "ak15NTracks", (200, 0, 200, 5)),

    ("SUEP_S1_SR", "boostedSphericity", (70, 0.3, 1.0, 5)),

    ("W_mt_SR", "wTransverseMass", (100, 30.0, 130.0, 5)),
    ("W_pt_SR", "wPt", (94, 60, 1000, 5)),
    ("W_phi_SR", "wPhi", (60, -3.2, 3.2, 3)),

    ("PuppiMET_pt_SR", "metPt", (94, 30, 500, 5)),
    ("PuppiMET_phi_SR", "metPhi", (100, -3, 3, 3)),
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
    return np.array(values[1:-1])

def get_hist_from_root(root_file, hist_name):
    for k in root_file.keys():
        if hist_name == k.split(";")[0]:
            return root_file[k].to_hist()
    return None

def scale_root(root_frac, gensumweight, root_norm_factor):
    root_raw_counts = root_frac * gensumweight
    root_vals = root_frac * root_norm_factor * gensumweight / gensumweight
    root_uncertainties = np.sqrt(root_raw_counts) / gensumweight * root_norm_factor
    return root_vals, root_uncertainties

def scale_saf(saf_data, norm_factor):
    saf_yields = saf_data * norm_factor
    saf_unc = np.sqrt(saf_data) * norm_factor
    return saf_yields, saf_unc

# 1D plotting function
def plot_1d_histos(saf_content, root_file, outdir, paired_histos, norm_factor, root_norm_factor, gensumweight):

    os.makedirs(outdir, exist_ok=True)

    for root_name, saf_name, (nbins, xmin, xmax, rebin) in paired_histos:

        saf_data = readSAFFile(saf_content, saf_name)
        if saf_data.size == 0:
            print(f"Missing SAF data for {saf_name}")
            continue

        root_hist = get_hist_from_root(root_file, root_name)
        if root_hist is None:
            print(f"Missing ROOT histogram {root_name}")
            continue

        # Original SAF bins and raw data
        orig_edges = np.linspace(xmin, xmax, nbins + 1)
        orig_centers = 0.5 * (orig_edges[:-1] + orig_edges[1:])
        # Rebin factor
        coarse_edges = orig_edges[::rebin]
        if coarse_edges[-1] != orig_edges[-1]:
            coarse_edges = np.append(coarse_edges, orig_edges[-1])
        bin_centers = 0.5 * (coarse_edges[:-1] + coarse_edges[1:])

        # Rebin SAF raw counts
        saf_raw = saf_data
        saf_coarse_raw = np.array([saf_raw[i:i+rebin].sum() for i in range(0, len(saf_raw), rebin)])
        saf_yields, saf_uncertainties = scale_saf(saf_coarse_raw, norm_factor)

        # Rebin ROOT raw counts via histogram
        root_frac = root_hist.values()
        root_edges = root_hist.axes[0].edges
        root_centers = 0.5 * (root_edges[:-1] + root_edges[1:])
        root_raw_counts = root_frac * gensumweight
        rebinned_raw, _ = np.histogram(root_centers, bins=coarse_edges, weights=root_raw_counts)
        root_vals = rebinned_raw / gensumweight * root_norm_factor
        root_uncertainties = np.sqrt(rebinned_raw) / gensumweight * root_norm_factor

        # Plot and ratio
        fig, (ax1, ax2) = plt.subplots(2, 1, gridspec_kw={'height_ratios': [3, 1]}, sharex=True, figsize=(7, 6))
        # Plot both on the rebinned bin centers
        ax1.errorbar(bin_centers, root_vals, yerr=root_uncertainties, fmt='o', color='darkgreen', label=f'CMSSW: {root_name}')
        ax1.errorbar(bin_centers, saf_yields, yerr=saf_uncertainties, fmt='o', color='navy', label=f'Delphes: {saf_name}')
        ax1.set_xlim(xmin, xmax)
        ax1.set_ylabel('Events')
        ax1.legend()
        ax1.set_title(f"{saf_name} vs {root_name}")

        # Ratio plot
        with np.errstate(divide='ignore', invalid='ignore'):
            ratio = np.true_divide(saf_yields, root_vals)
            ratio[~np.isfinite(ratio)] = 0
        ax2.step(bin_centers, ratio, where='mid', color='black')
        ax2.axhline(1, color='gray', linestyle='--')
        ax2.set_ylabel('Ratio')
        ax2.set_xlabel(saf_name)
        ax2.set_ylim(0, 2)

        plt.tight_layout()
        out_pdf = os.path.join(outdir, f"{saf_name}.pdf")
        out_png = os.path.join(outdir, f"{saf_name}.png")
        plt.savefig(out_pdf)
        plt.savefig(out_png)
        plt.close()
        print(f"Saved 1D plot: {out_png}")

# 2D slicing function: compares 2D ROOT histogram slices to SAF ABCD regions
def plot_2d_slices(saf_content, root_file, outdir, norm_factor, root_norm_factor, gensumweight):
    # ABCD regions and y‐ranges
    region_names = ['ABCD_SR0','ABCD_SR1','ABCD_SR2','ABCD_SR3','ABCD_SR4']
    y_ranges     = [(30,40),(40,50),(50,60),(60,80),(80,120)]

    # Access 2D ROOT histogram
    h2 = get_hist_from_root(root_file, '2D_SUEP_S1_vs_SUEP_nconst_SR')
    if h2 is None:
        raise KeyError('2D_SUEP_S1_vs_SUEP_nconst_SR not found in ROOT file')
    x_edges = h2.axes[0].edges
    y_edges = h2.axes[1].edges
    vals2d  = h2.values()

    # Compute centers and mask for x in [0.5,1.0)
    x_centers = 0.5*(x_edges[:-1]+x_edges[1:])
    y_centers = 0.5*(y_edges[:-1]+y_edges[1:])
    x_mask    = (x_centers>=0.5)&(x_centers<1.0)

    slice_centers = []
    root_counts   = []
    root_unc      = []
    saf_counts    = []
    saf_unc       = []

    for name,(y_low,y_high) in zip(region_names,y_ranges):
        y_mask = (y_centers>=y_low)&(y_centers<y_high)
        mask2d = np.outer(x_mask,y_mask)
        raw2d  = np.sum(vals2d[mask2d]*gensumweight)
        scaled = raw2d/gensumweight*root_norm_factor
        unc    = np.sqrt(raw2d)/gensumweight*root_norm_factor
        root_counts.append(scaled)
        root_unc.append(unc)
        # SAF summation using readSAFFile
        saf_vals = readSAFFile(saf_content,name)
        raw1d    = np.sum(saf_vals)
        scaled1d = raw1d*norm_factor
        unc1d    = np.sqrt(raw1d)*norm_factor
        saf_counts.append(scaled1d)
        saf_unc.append(unc1d)
        slice_centers.append((y_low+y_high)/2)

    # Plot
    fig,ax = plt.subplots()
    ax.errorbar(slice_centers,root_counts,yerr=root_unc,fmt='o',color='darkgreen',label='CMSSW 2D')
    ax.errorbar(slice_centers,saf_counts, yerr=saf_unc, fmt='o',color='navy',    label='Delphes ABCD')
    ax.set_xticks(slice_centers)
    ax.set_xticklabels(region_names)
    ax.set_xlabel('ABCD Region')
    ax.set_ylabel('Scaled yield')
    ax.legend()
    plt.tight_layout()
    out_png = os.path.join(outdir,'ABCD_SR.png')
    plt.savefig(out_png)
    plt.close()
    print(f"Saved 2D slices plot: {out_png}")

def main():
    parser = argparse.ArgumentParser(description="Overlay SAF vs ROOT histograms with ratio")
    parser.add_argument("--saf", required=True, help="SAF file")
    parser.add_argument("--root", required=True, help="ROOT file")
    parser.add_argument("--outdir", default="ratio_plots", help="Output directory")
    parser.add_argument("--total_gen", type=float, required=True, help="Total generated events")
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    with open(args.saf, "r") as f:
        saf_content = f.read()

    root_file = uproot.open(args.root)
    tfile_py = rpy.TFile.Open(args.root)     # Read metadata using PyROOT because uproot's TObjString model doesn't expose the string member

    mddir = tfile_py.Get("metadata")
    if not mddir:
        raise KeyError("Metadata directory not found in ROOT file via PyROOT either")
    
    # Extract string values and convert to float (convert TString to Python str first)
    gensumweight = float(str(mddir.Get("gensumweight").GetString()))
    tfile_py.Close()

    xsec = 0.09426 * 3 * 1373  # (WH process σ×BR; see https://twiki.cern.ch/twiki/bin/view/LHCPhysics/CERNYellowReportPageAt13TeV#WH_Process)
    lumi = 59.9 # 2018 Lumi

    # SAF normalization
    norm_factor = lumi * xsec / args.total_gen

    # ROOT normalization
    root_norm_factor = lumi * xsec

    # generate plots
    plot_1d_histos(saf_content, root_file, args.outdir, paired_histos, norm_factor, root_norm_factor, gensumweight)
    plot_2d_slices(saf_content, root_file, args.outdir, norm_factor, root_norm_factor, gensumweight)

if __name__ == "__main__":
    main()

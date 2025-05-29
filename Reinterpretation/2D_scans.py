import argparse
import os
import numpy as np
import matplotlib.pyplot as plt

# Axis limits
binning_dict = {
    "lepPt": (30, 400),
    "lepEta": (-2.5, 2.5),
    "ak15NTracks": (0.0, 75.0),
    "wTransverseMass": (30.0, 130.0),
    "wPt": (60, 1000),
    "ak4Pt": (0, 1000),
    "boostedSphericity_CENTRAL": (0, 1.0),
    "metPt": (30, 500),
    "NJets": (0, 20),
    "ak15Pt": (0, 300),
    "wPhi": (-3.2, 3.2),
    "W_SUEP_ratio": (0, 3)
}

def manual_pearsonr(x, y):
    x_mean = np.mean(x)
    y_mean = np.mean(y)
    cov = np.sum((x - x_mean) * (y - y_mean))
    std_x = np.sqrt(np.sum((x - x_mean) ** 2))
    std_y = np.sqrt(np.sum((y - y_mean) ** 2))
    return cov / (std_x * std_y)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--csv", required=True)
    parser.add_argument("--outdir", default="plots")
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    # Load data from CSV
    data = np.genfromtxt(args.csv, delimiter=",", names=True)

    # Print event numbers and ak15Pt for ak15Pt < 60
    if "ak15Pt" in data.dtype.names and "EventNumber" in data.dtype.names:
        low_pt_mask = data["ak15Pt"] < 60
        event_numbers = data["EventNumber"][low_pt_mask]
        ak15pts = data["ak15Pt"][low_pt_mask]

        print("Events with ak15Pt < 60:")
        for evt, pt in zip(event_numbers, ak15pts):
            print(f"Event {int(evt)}: ak15Pt = {pt:.2f}")
            print(f"\nTotal events with ak15Pt < 60: {len(event_numbers)}")
    else:
        print("Missing 'ak15Pt' or 'EventNumber' in CSV columns.")

    yname = "ak15NTracks"
    y_vals = data[yname]

    for xname in binning_dict:
        if xname == yname or xname not in data.dtype.names:
            continue

        x_vals = data[xname]

        # Filter finite values
        mask = np.isfinite(x_vals) & np.isfinite(y_vals)
        x_vals_clean = x_vals[mask]
        y_vals_clean = y_vals[mask]

        # Pearson correlation (manual)
        corr = manual_pearsonr(x_vals_clean, y_vals_clean)

        # Axis binning
        x_min, x_max = binning_dict[xname]
        y_min, y_max = binning_dict[yname]
        x_bins = 50
        y_bins = 50

        # Plot
        plt.figure(figsize=(8, 6))
        plt.hist2d(x_vals_clean, y_vals_clean,
                   bins=[np.linspace(x_min, x_max, x_bins + 1),
                         np.linspace(y_min, y_max, y_bins + 1)],
                   cmap='viridis')
        plt.colorbar(label="Counts")
        plt.xlabel(xname)
        plt.ylabel(yname)
        plt.title(f"{xname} vs {yname}")
        plt.text(0.05, 0.95, f"Pearson r = {corr:.3f}",
                 transform=plt.gca().transAxes,
                 fontsize=12, verticalalignment='top',
                 bbox=dict(boxstyle="round", facecolor="white", alpha=0.8))
        plt.tight_layout()
        outpath = os.path.join(args.outdir, f"{xname}_vs_{yname}.png")
        plt.savefig(outpath)
        plt.close()
        print(f"Saved: {outpath}")

if __name__ == "__main__":
    main()

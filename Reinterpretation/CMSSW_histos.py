import sys
import uproot
import matplotlib.pyplot as plt
import os
import gc

infile_name = sys.argv[1]
output_dir = sys.argv[2]

# Make output directory if it doesn't exist
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

# Load ROOT histograms
# source code: https://github.com/SUEPPhysics/SUEPCoffea_dask/blob/342406f140ae4a1cb8d9c86ab3bcb8abc1be1551/plotting/utils/loader.py#L426
def openroot(infile_name):
    _plots = {}
    _metadata = {}
    with uproot.open(infile_name) as _infile:
        for k in _infile.keys():
            if "metadata" == k.split(";")[0]:
                for kk in _infile[k].keys():
                    _metadata[kk.split(";")[0]] = _infile[k][kk].title()
            elif "metadata" not in k:
                _plots[k.split(";")[0]] = _infile[k].to_hist()
    gc.collect()
    return _plots, _metadata

# Load histograms
plots, metadata = openroot(infile_name)

# Plot and save each histogram
for name, hist in plots.items():
    fig, ax = plt.subplots()
    hist.plot(ax=ax)

    # Add metadata to plot title if available
    title = metadata.get('title', name)
    ax.set_title(title)

    output_file = os.path.join(output_dir, f"{name}.pdf")
    plt.savefig(output_file)
    plt.close(fig)

    print(f"Saved plot: {output_file}")

print("Finished plotting all histograms.")

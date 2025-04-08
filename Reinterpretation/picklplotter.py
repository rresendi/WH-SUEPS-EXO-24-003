import pickle
# import matplotlib.pyplot as plt

# Function to open and read pickle file
def openpickle(infile_name):
    _plots = {}
    _metadata = {}
    with open(infile_name, "rb") as openfile:
        while True:
            try:
                input = pickle.load(openfile)
                _plots.update(input["hists"].copy())
                _metadata.update(input["metadata"].copy())
            except EOFError:
                break
    return _plots, _metadata

# Load data from pickle file
plots, metadata = openpickle('CMSSW_HISTOS/forRyleigh/WH_3_6_signal_2018/SUEP_mS125.000_mPhi3.000_T3.000_modegeneric.pkl')
print(plots)
# Iterate over histograms and plot them
for hist_name, (bin_counts, bin_edges) in plots.items():
    plt.figure(figsize=(8, 6))

    # set bin edges
    plt.hist(bin_edges[:-1], bins=bin_edges, weights=bin_counts, 
             histtype='stepfilled', label=hist_name, color = 'navy')

    # add metadata
    plt.title(metadata.get('title', f'Histogram: {hist_name}'))
    plt.xlabel(f'{hist_name}')
    plt.ylabel(metadata.get('ylabel', 'Events'))

    plt.legend()
    plt.tight_layout()

    # Display the plot
    outname = f"CMSSW{hist_name}.pdf"
    plt.savefig()

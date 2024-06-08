import uproot
import os
import argparse
import awkward as ak
import numpy as np
import ROOT
from array import array

# Sets batch mode so no popup window                                                                                                                                                                       
ROOT.gROOT.SetBatch(True)

# Initialize parser                                                                                                                                                                                        
parser = argparse.ArgumentParser()
parser.add_argument("--input", help="Name of input file", type=str)
args = vars(parser.parse_args())

# Name of sample                                                                                                                                                                                           
sample_name = args["input"]

output_file = "MC_electron_efficiencies.root"
input_file = "/eos/user/j/jreicher/SUEP/WH_private_signals/merged/" + sample_name + ".root"

# suep decay type                                                                                                                                                                                          
if "generic" in sample_name:
    decay_type = "generic"
elif "hadronic" in sample_name:
    decay_type = "hadronic"
else:
    decay_type = "leptonic"

# conditions for what year                                                                                                                                                                                 
if "UL18" in sample_name:
    year = "2018 conditions"    
    folder = "ele_eff_outputs_2017/"
elif "UL16APV" in sample_name:
    year = "2016 APV conditions"
    folder = "ele_eff_outputs_2016APV/"
else:
    year = "2016 conditions"
    folder = "ele_eff_outputs_2016/"

# dark meson (phi) mass                                                                                                                                                                                    
if "MD2.00" in sample_name:
    md = "2.00 [GeV]"
elif "MD4.00" in sample_name:
    md = "4.00 [GeV]"
elif "MD3.00" in sample_name:
    md = "3.00 [GeV]"
elif "MD8.00" in sample_name:
    md = "8.00 [GeV]"
elif "MD1.00" in sample_name:
    md = "1.00 [GeV]"
else:
    md = "1.40 [GeV]"

# temperature                                                                                                                                                                                              
if "T0.25" in sample_name:
    temp = "0.25"
if "T0.35" in sample_name:
    temp = "0.35"
if "T0.50" in sample_name:
    temp = "0.50"
elif "T0.75" in sample_name:
    temp = "0.75"
elif "T1.00" in sample_name:
    temp = "1.00"
elif "T1.50" in sample_name:
    temp = "1.50"
elif "T2.00" in sample_name:
    temp = "2.00"
elif "T3.00" in sample_name:
    temp = "3.00"
elif "T4.00" in sample_name:
    temp = "4.00"
elif "T8.00" in sample_name:
    temp = "8.00"
elif "T12.00" in sample_name:
    temp = "12.00"
elif "T16.00" in sample_name:
    temp = "16.00"
elif "T32.00" in sample_name:
    temp = "32.00"
else:
    temp = "6.00"

# Gets relevant events from file                                                                                                                                                                           
def Events(f):
    evs = f['Events'].arrays(['HLT_Ele32_WPTight_Gsf',
                              'HLT_Ele115_CaloIdVT_GsfTrkIdT',
                              'HLT_Photon175',
                              'HLT_Photon200',
                              'Electron_cutBased',
                              'Electron_pt',
                              'Electron_mvaFall17V2Iso_WP80',
                              'Electron_eta',
                              'Electron_dxy',
                              'Electron_dz',
                              'Electron_phi',
                              'Electron_mass',
                              'Electron_pdgId',
                              'TrigObj_pt',
                              'TrigObj_eta',
                              'TrigObj_phi',
                              'TrigObj_id',
                              'TrigObj_filterBits'
                              ])
    return evs

with uproot.open(input_file) as f:
    evs = Events(f)

# Defining a good electron                                                                                                                                                                                \                                                                                                                                                                                       

electrons = ak.zip({
        "pt": evs["Electron_pt"],
        "eta": evs["Electron_eta"],
        "phi": evs["Electron_phi"],
        "mass": evs["Electron_mass"],
        "charge": evs["Electron_pdgId"]/(-11),
        "pdgId": evs["Electron_pdgId"],
        "isTight": evs["Electron_mvaFall17V2Iso_WP80"],
        "isTightIso": evs["Electron_mvaFall17V2Iso_WP80"]
        }, with_name = "Momentum4D")

cutElectrons = (
        (evs["Electron_cutBased"] >= 2)
        & (evs["Electron_pt"] >= 15)
        & (evs["Electron_mvaFall17V2Iso_WP80"])
        & (abs(evs["Electron_dxy"]) < 0.05 + 0.05 * (abs(evs["Electron_eta"]) > 1.479))
        & (abs(evs["Electron_dz"]) < 0.10 + 0.10 * (abs(evs["Electron_eta"]) > 1.479))
        & ((abs(evs["Electron_eta"]) < 1.444) | (abs(evs["Electron_eta"]) > 1.566))
        & (abs(evs["Electron_eta"]) < 2.5)
    )

goodElectrons = electrons[cutElectrons]

# Computes deltaR2
def deltaR2(eta1, phi1, eta2, phi2):
    deta = eta1 - eta2
    dphi = np.arccos(np.cos(phi1 - phi2))

    return deta ** 2 + dphi ** 2

# Setting up alldr2 to be usable

alldr2 = None

# Gets matched online/offline electrons from file
def isHLTMatched(events, electrons):
    
    global alldr2
    
    trigObj = ak.zip({
        "pt": events["TrigObj_pt"],
        "eta": events["TrigObj_eta"],
        "phi": events['TrigObj_phi'],
        "mass": 0.,
        "id": events['TrigObj_id'],
        "filterBits": events['TrigObj_filterBits']
        }, with_name="Momentum4D")

    filterbits1 = ((events['TrigObj_filterBits'] & 2) == 2)
    filterbits2 = ((events['TrigObj_filterBits'] & 2048) == 2048)
    filterbits3 = ((events['TrigObj_filterBits'] & 8192) == 8192)
    
    trigObjSingleEl = trigObj[((abs(trigObj.id) == 11)
                               & (trigObj.pt >= 35)
                               & ((abs(trigObj.eta) < 1.444) | (abs(trigObj.eta) > 1.566))
                               & (abs(trigObj.eta) < 2.5)
                               & (filterbits1 |
                                filterbits2 |
                                filterbits3))]
    
    toMatch1El, trigObjSingleEl = ak.unzip(ak.cartesian([goodElectrons, trigObjSingleEl], axis=1, nested=True))
    alldr2 = deltaR2(toMatch1El.eta, toMatch1El.phi, trigObjSingleEl.eta, trigObjSingleEl.phi)
    match1El = (ak.sum(ak.where(ak.min(alldr2, axis=2) < 0.1, True, False), axis=1) >= 1)

    return match1El

# Defines binning and histograms
ele_bin_edges = array('d', [0, 2, 4, 6, 8, 10, 12,
                            14, 16, 18, 20, 22,
                            24, 26, 28, 30, 32,
                            34, 36, 38, 40, 50,
                            60, 70, 80, 90, 100,
                            120, 140, 160, 180, 200])

# Create a 2D histogram for deltaR vs. pT for each eta bin                                                                                                                                                 
deltaR_vs_pt_eta1 = ROOT.TH1D("deltaR_vs_pt_eta1", "DeltaR vs pT (|eta| < 1)",len(ele_bin_edges)-1,ele_bin_edges)
deltaR_vs_pt_eta2 = ROOT.TH1D("deltaR_vs_pt_eta2", "DeltaR vs pT (1 <= |eta| < 2)",len(ele_bin_edges)-1,ele_bin_edges)
deltaR_vs_pt_eta3 = ROOT.TH1D("deltaR_vs_pt_eta3", "DeltaR vs pT (2 <= |eta| < 2.5)",len(ele_bin_edges)-1,ele_bin_edges)

# Function for filling the deltaR histograms and 2D histograms                                                                                                                                             
def ele_hists(events, etas, hists):
    eta_min = etas[0]
    eta_max = etas[1]

    # Electron selection                                                                                                                                                                                   
    ele_quality_check = isHLTMatched(events, goodElectrons)

    # Trigger selection                                                                                                                                                                                    

    if "UL17" in sample_name or "UL18" in sample_name:
        triggerSingleElectron = (
            events["HLT_Ele32_WPTight_Gsf"] |
            events["HLT_Ele115_CaloIdVT_GsfTrkIdT"] |
            events["HLT_Photon200"]
        )
    else:
        triggerSingleElectron = (
            events["HLT_Ele32_WPTight_Gsf"] |
            events["HLT_Ele115_CaloIdVT_GsfTrkIdT"] |
            events["HLT_Photon175"]
        )
    
    # Cut on eta                                                                                                                                                                                           
    eta_split = (
        (np.abs(events["Electron_eta"]) >= eta_min) &
        (np.abs(events["Electron_eta"]) < eta_max)
    )
    
    
    # Compute deltaR and fill 2D histograms, basically i need to extraoplate alldr2 and plot that                                                                                                                                                              
    trigObj = ak.zip({
        "eta": events["TrigObj_eta"],
        "phi": events['TrigObj_phi']
    }, with_name = "Momentum4D")
    trigObjSingleEl = trigObj[((abs(events['TrigObj_id']) == 11) &
                               ((events['TrigObj_filterBits'] & 2) |
                                (events['TrigObj_filterBits'] & 2048) |
                                (events['TrigObj_filterBits'] & 8192)))]
    toMatch1El, trigObjSingleEl = ak.unzip(ak.cartesian([goodElectrons, trigObjSingleEl], axis = 1, nested = True))
    alldr2 = deltaR2(toMatch1El.eta, toMatch1El.phi, trigObjSingleEl.eta, trigObjSingleEl.phi)
    deltaR_values = ak.min(alldr2, axis=2)
    
    
    ele = events["Electron_pt"]
    evs = ele[ele_quality_check & eta_split]
    tr_evs = evs[triggerSingleElectron]

    for deltaR, pt in zip(deltaR_values[both], ele[both]):
        for value, pt_val in zip(deltaR, pt):
            deltaR_vs_pt_hist.Fill(pt_val, value)

# Fill 2D histograms                                                                                                                                                                                       
eta_split = [[0, 1.0], [1.0, 2.0], [2.0, 3.0]]
eta_hists = [[deltaR_vs_pt_eta1],
             [deltaR_vs_pt_eta2],
             [deltaR_vs_pt_eta3]]


for (etas, hists) in zip(eta_split, eta_hists):
    ele_hists(evs, etas, hists)

# Create canvas for plotting                                                                                                                                                                               
c1 = ROOT.TCanvas("canvas", "", 800, 600)

# Creates deltaR Plot w legend                                                                                                                                                                         

eta1_effs.SetTitle("deltaR values of cut events in bins of pT;Electron pT [GeV];deltaR")
legend = ROOT.TLegend(0.5,0.1,0.9,0.4)
legend.AddEntry(deltaR_vs_pt_eta1, "|#eta|<1.0", "l")
legend.AddEntry(deltaR_vs_pt_eta2, "1.0<|#eta|<2.0", "l")
legend.AddEntry(deltaR_vs_pt_eta3, "2.0<|#eta|<3.0", "l")
legend.AddEntry(ROOT.nullptr, "T = " + temp + "GeV, " + year,"")
legend.AddEntry(ROOT.nullptr, "SUEP decay type: " + decay_type,"")
legend.AddEntry(ROOT.nullptr, "Dark meson mass = " + md,"")
legend.SetTextColor(ROOT.kBlack)
legend.SetTextFont(42)
legend.SetTextSize(0.03)

# Draw plot                                                                                                                                                                                                

eta1_effs.Draw()
eta2_effs.SetLineColor(ROOT.kRed)
eta2_effs.Draw("same")
eta3_effs.SetLineColor(ROOT.kBlue)
eta3_effs.Draw("same")
legend.Draw("same")
c1.Update()

# Saves to pdf                                                                                                                                                                                             

c1.SaveAs(folder + sample_name + "_deltaR.pdf")

# Saves overall efficiency

root_file = ROOT.TFile(output_file, "UPDATE")
root_file.cd()

eff_dir = root_file.Get("Efficiencies")
if not eff_dir:
    eff_dir = root_file.mkdir("Efficiencies")
eff_dir.cd()
ele_eff.Write()

root_file.Close()

print("sample " + sample_name + " complete")

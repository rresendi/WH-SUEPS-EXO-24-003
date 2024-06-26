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

parser.add_argument("--input", help = "Name of input file", type = str)
args = vars(parser.parse_args())

# Name of sample                                                                                                                                                                                           
sample_name = args["input"]

output_file = "MC_electron_efficiencies.root"
input_file = "/eos/user/j/jreicher/SUEP/WH_private_signals/merged/" + sample_name + ".root"

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
                ])
    return evs

with uproot.open(input_file) as f:
    evs = Events(f)

# Defines binning and histograms                                                                                                                                                                           

ele_bin_edges = array('d',[0,2,4,6,8,10,12,
                         14,16,18,20,22,
                         24,26,28,30,32,
                         34,36,38,40,50,
                         60,70,80,90,100,
                         120,140,160,180,200])

# Histograms for overall efficiency                                                                                                                                                                        

ele_totalhist = ROOT.TH1D("total_events","Total Events",len(ele_bin_edges)-1,ele_bin_edges)
ele_filthist = ROOT.TH1D("filt_events","Filtered Events",len(ele_bin_edges)-1,ele_bin_edges)

# Split into three regions of eta                                                                                                                                                                          

eta1_ele_totalhist = ROOT.TH1D("total_events","Total Events",len(ele_bin_edges)-1,ele_bin_edges)
eta1_ele_filthist = ROOT.TH1D("filt_events","Filtered Events",len(ele_bin_edges)-1,ele_bin_edges)
eta2_ele_totalhist = ROOT.TH1D("total_events","Total Events",len(ele_bin_edges)-1,ele_bin_edges)
eta2_ele_filthist = ROOT.TH1D("filt_events","Filtered Events",len(ele_bin_edges)-1,ele_bin_edges)
eta3_ele_totalhist = ROOT.TH1D("total_events","Total Events",len(ele_bin_edges)-1,ele_bin_edges)
eta3_ele_filthist = ROOT.TH1D("filt_events","Filtered Events",len(ele_bin_edges)-1,ele_bin_edges)
eta4_ele_totalhist = ROOT.TH1D("total_events","Total Events",len(ele_bin_edges)-1,ele_bin_edges)
eta4_ele_filthist = ROOT.TH1D("filt_events","Filtered Events",len(ele_bin_edges)-1,ele_bin_edges)

# Function for filling the histograms                                                                                                                                                                      

def ele_hists(events, etas, hists):
    ele_totalhist = hists[0]
    ele_filthist = hists[1]
    eta_min = etas[0]
    eta_max = etas[1]

    # Electron selection                                                                                                                                                                                   

    ele_quality_check = (
        (evs["Electron_cutBased"] >= 2)
        & (evs["Electron_mvaFall17V2Iso_WP80"])
        & (abs(evs["Electron_dxy"]) < 0.05 + 0.05 * (abs(evs["Electron_eta"]) > 1.479))
        & (abs(evs["Electron_dz"]) < 0.10 + 0.10 * (abs(evs["Electron_eta"]) > 1.479))
        & ((abs(evs["Electron_eta"]) < 1.444) | (abs(evs["Electron_eta"]) > 1.566))
        & (abs(evs["Electron_eta"]) < 2.5)
    )
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
        (np.abs(events["Electron_eta"]) >= eta_min)
        & (np.abs(events["Electron_eta"]) < eta_max)
    )

    # Select based on trigger                                                                                                                                                                              

    ele = events["Electron_pt"]
    evs = ele[ele_quality_check & eta_split]
    tr_evs = evs[triggerSingleElectron]

    # Fill histograms                                                                                                                                                                                      

    for ev in evs:
        for entry in ev:
            ele_totalhist.Fill(entry)
    for ev in tr_evs:
        for entry in ev:
            ele_filthist.Fill(entry)

    return 0

eta_split = [[0, 1.0], [1.0, 2.0], [2.0, 3.0]]
eta_hists = [[eta1_ele_totalhist,eta1_ele_filthist],
             [eta2_ele_totalhist,eta2_ele_filthist],
             [eta3_ele_totalhist,eta3_ele_filthist]]

for (etas,hists) in zip(eta_split, eta_hists):
    ele_hists(evs, etas, hists)

# Fills efficiency                                                                                                                                                                                         

eta1_effs = ROOT.TEfficiency(eta1_ele_filthist,eta1_ele_totalhist)
eta2_effs = ROOT.TEfficiency(eta2_ele_filthist,eta2_ele_totalhist)
eta3_effs = ROOT.TEfficiency(eta3_ele_filthist,eta3_ele_totalhist)
c1 = ROOT.TCanvas ("canvas","",800,600)

# Get overall Efficiency:                                                                                                                                                                                  

ele_eff = ele_filthist.Clone()

ele_eff.SetName(sample_name)
ele_eff.Sumw2()
ele_eff.Divide(ele_totalhist)

# Creates Efficiency Plot w legend                                                                                                                                                                         

eta1_effs.SetTitle("Electron Trigger Efficiency in bins of pT;Electron pT [GeV];Efficiency")
legend = ROOT.TLegend(0.5,0.1,0.9,0.4)
legend.AddEntry(eta1_effs, "|#eta|<1.0", "l")
legend.AddEntry(eta2_effs, "1.0<|#eta|<2.0", "l")
legend.AddEntry(eta3_effs, "2.0<|#eta|<3.0", "l")
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

c1.SaveAs(folder + sample_name + "_Efficiency.pdf")

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

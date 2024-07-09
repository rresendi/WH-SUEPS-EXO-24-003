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
    folder = "ele_eff_outputs_2018/"
elif "UL17" in sample_name:
    year = "2017 conditions"
    folder = "ele_eff_outputs_2017/"
elif "UL16APV" in sample_name:
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

# Function for filling the histograms                                                                                                                                                                      

def ele_hists(events, etas, hists):
    ele_totalhist = hists[0]
    ele_filthist = hists[1]
    eta_min = etas[0]
    eta_max = etas[1]
    
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
    
    # quality requirements for electrons                                                                                                                                                                              
    ele_quality_check = (
            (events["Electron_cutBased"] >= 2)
            & (events["Electron_mvaFall17V2Iso_WP80"])
            & (abs(events["Electron_dxy"]) < 0.05 + 0.05 * (abs(events["Electron_eta"]) > 1.479))
            & (abs(events["Electron_dz"]) < 0.10 + 0.10 * (abs(events["Electron_eta"]) > 1.479))
            & ((abs(events["Electron_eta"]) < 1.444) | (abs(events["Electron_eta"]) > 1.566))
            & (abs(events["Electron_eta"]) < 2.5)
            )
    
    # Cut on eta                                                                                                                                                                                           

    eta_split = (
        (np.abs(events["Electron_eta"]) >= eta_min)
        & (np.abs(events["Electron_eta"]) < eta_max)
    )
    
    # Combined mask
    
    combined_mask = ele_quality_check & eta_split
    
    # Select based on trigger                                                                                                                                                                              

    ele = events["Electron_pt"]
    evs = ele[eta_split]
    tr_evs = evs[combined_mask & triggerSingleElectron]

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

# Loop over eta bins and fill histograms for numerator and denominator plots
ele_hists(evs, eta_split[0], [eta1_ele_totalhist, et1a_ele_filthist])
ele_hists(evs, eta_split[1], [eta2_ele_totalhist, eta2_ele_filthist])
ele_hists(evs, eta_split[2], [eta3_ele_totalhist, eta3_ele_filthist])

# Fills efficiency                                                                                                                                                                                         

eta1_effs = ROOT.TEfficiency(eta1_ele_filthist,eta1_ele_totalhist)
eta2_effs = ROOT.TEfficiency(eta2_ele_filthist,eta2_ele_totalhist)
eta3_effs = ROOT.TEfficiency(eta3_ele_filthist,eta3_ele_totalhist)
c1 = ROOT.TCanvas ("canvas","",800,600)

# Get overall Efficiency:                                                                                                                                                                                  

ele_eff = ele_filthist.Clone()

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
eta1_effs.Draw("AP")

# Update canvas                                                                                                                                                                           
c1.Update()

# Get painted graph and set y-axis range                                                                                                                                                  
efficiency = eta1_effs
graph = efficiency.GetPaintedGraph()
graph.SetMinimum(0)

# Update canvas again                                                                                                                                                                     
ROOT.gPad.Update()

# Draw legend                                                                                                                                                                             
legend.Draw()
c1.Update()
efficiency = eta1_effs
efficiency.Draw()
ROOT.gPad.Update()
graph = efficiency.GetPaintedGraph()
graph.SetMinimum(0)

ROOT.gPad.Update()
legend.Draw()
c1.Update()
eta2_effs.Draw("P same")
eta2_effs.SetLineColor(ROOT.kRed)
eta3_effs.Draw("P same")
eta3_effs.SetLineColor(ROOT.kBlue)

# Saves to pdf

c1.SaveAs(sample_name + "_Efficiency.pdf")

# Set line colors for histograms                                                                                                                                                                           
eta1_ele_filthist.SetLineColor(ROOT.kBlack)
eta2_ele_filthist.SetLineColor(ROOT.kRed)
eta3_ele_filthist.SetLineColor(ROOT.kBlue)
eta1_ele_totalhist.SetLineColor(ROOT.kBlack)
eta2_ele_totalhist.SetLineColor(ROOT.kRed)
eta3_ele_totalhist.SetLineColor(ROOT.kBlue)

# Create canvas for numerator plot                                                                                                                                                                         
c2 = ROOT.TCanvas("canvas1", "", 800, 600)

# Create first legend                                                                                                                                                                                      
legend = ROOT.TLegend(0.5, 0.6, 0.9, 0.9)
legend.AddEntry(eta1_ele_filthist, "|#eta|<1.0", "l")
legend.AddEntry(eta2_ele_filthist, "1.0<|#eta|<2.0", "l")
legend.AddEntry(eta3_ele_filthist, "2.0<|#eta|<3.0", "l")
legend.AddEntry(ROOT.nullptr, "T = " + temp + "GeV, " + year,"")
legend.AddEntry(ROOT.nullptr, "SUEP decay type: " + decay_type,"")
legend.AddEntry(ROOT.nullptr, "Dark meson mass = " + md,"")
legend.SetTextColor(ROOT.kBlack)
legend.SetTextFont(42)
legend.SetTextSize(0.03)

# Plot numerator for eta bins 1, 2, and 3                                                                                                                                                                  
c2.cd(1)
eta1_ele_filthist.SetTitle("Number of Cut Events in bins of pT;Electron pT [GeV];Number of events")
eta1_ele_filthist.SetMaximum(550)
eta1_ele_filthist.SetStats(0)
eta1_ele_filthist.Draw()
eta2_ele_filthist.Draw("same")
eta3_ele_filthist.Draw("same")
legend.Draw()

# Save numerator plot as PDF                                                                                                                                                                               
c2.SaveAs(sample_name + "_numerator_eta_bins.pdf")

# Create canvas for denominator plot                                                                                                                                                                       
c3 = ROOT.TCanvas("canvas2", "", 800, 600)

# Create second legend                                                                                                                                                                                     
legend = ROOT.TLegend(0.5, 0.6, 0.9, 0.9)
legend.AddEntry(eta1_ele_filthist, "|#eta|<1.0", "l")
legend.AddEntry(eta2_ele_filthist, "1.0<|#eta|<2.0", "l")
legend.AddEntry(eta3_ele_filthist, "2.0<|#eta|<3.0", "l")
legend.AddEntry(ROOT.nullptr, "T = " + temp + "GeV, " + year,"")
legend.AddEntry(ROOT.nullptr, "SUEP decay type: " + decay_type,"")
legend.AddEntry(ROOT.nullptr, "Dark meson mass = " + md,"")
legend.SetTextColor(ROOT.kBlack)
legend.SetTextFont(42)
legend.SetTextSize(0.03)

# Plot denominator for eta bins 1, 2, and 3                                                                                                                                                                
c3.cd(1)
eta1_ele_totalhist.SetTitle("Number of Total Events in bins of pT;Electron pT [GeV];Number of Events")
eta1_ele_totalhist.SetStats(0)
eta1_ele_totalhist.Draw()
eta2_ele_totalhist.Draw("same")
eta3_ele_totalhist.Draw("same")
legend.Draw()
# Save denominator plot as PDF                                                                                                                                                                             
c3.SaveAs(sample_name + "_denominator_eta_bins.pdf")

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

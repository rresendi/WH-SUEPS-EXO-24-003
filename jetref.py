import ROOT
import os, sys
import argparse
import numpy as np
from array import array
import json

# Define lepton type
lepton = sys.argv[1] # Electron or Muon

# Data ?
data = sys.argv[2] # mc or data

# Eta option: Eta or noEta
etaOption = sys.argv[3]

# Data era
era = sys.argv[4] # 2016, 2016APV, 2017, 2018

lep1pt_bin_edges = array('d',[0, 2, 4, 6, 8, 10, 12,
                             14, 16, 18, 20, 22,
                            24, 26, 28, 30, 32,
                           34, 36, 38, 40, 50,
                          60, 70, 80, 90, 100,
                         120, 140, 160, 180, 200])

refhlt = "HLT_PFMETNoMu120_PFMHTNoMu120_IDTight"

# Lepton-specific configurations
if lepton == "Muon":
    inputFiles = sys.argv[6:]
    outputHistos = sys.argv[5] # "muon_efficiencies.root"
    hlt = ["HLT_IsoMu27", "HLT_Mu50"]
    offlineCuts = {
        "lep1pt": 40,
        "MET": 40,
        "mT" : (30, 130)
    }

    # Update histBins dictionary
    histBins = {
        "lep1pt": lep1pt_bin_edges,
        "MET": [30, 0, 300],
        "mT": [15, 0, 150],
        "lep1phi": [35, 0, 3.5],
        "MET phi": [35, 0, 3.5]
    }

    pt_label = "Muon pT [GeV]"
    eta_bins = ["eta1", "eta2", "eta3"]
    eta_ranges = [(0.0, 0.9), (0.9, 2.1), (2.1, 2.4)]
else:  # Electron-specific configurations
    inputFiles = sys.argv[6:]
    outputHistos = sys.argv[5] # "electron_efficiencies.root"
    if era == "2017" or era == "2018":
        hlt = ["HLT_Ele32_WPTight_Gsf", "HLT_Ele115_CaloIdVT_GsfTrkIdT", "HLT_Photon200"]
    else:
        hlt = ["HLT_Ele32_WPTight_Gsf", "HLT_Ele115_CaloIdVT_GsfTrkIdT", "HLT_Photon175"]
    offlineCuts = {
        "lep1pt": 30,
        "MET": 40,
        "mT" : (30, 130)
    }

    # Update histBins dictionary
    histBins = {
        "lep1pt": lep1pt_bin_edges,
        "MET": [30, 0, 300],
        "mT": [15, 0, 150],
        "lep1phi": [35, 0, 3.5],
        "MET phi": [35, 0, 3.5]
    }

    pt_label = "Electron pT [GeV]"
    eta_bins = ["eta1", "eta2", "eta3"]
    eta_ranges = [(0.0, 1.0), (1.0, 2.0), (2.0, 3.0)]

# Define reference cut function
def passRefCut(ev, era):
    if era in ["2018", "2017"]:
        cutAnyFilter = (
            ev.Flag_goodVertices
            and ev.Flag_globalSuperTightHalo2016Filter
            and ev.Flag_HBHENoiseFilter
            and ev.Flag_HBHENoiseIsoFilter
            and ev.Flag_EcalDeadCellTriggerPrimitiveFilter
            and ev.Flag_BadPFMuonFilter
            and ev.Flag_BadPFMuonDzFilter
            and ev.Flag_eeBadScFilter
            and ev.Flag_ecalBadCalibFilter
        )
    elif era in ["2016", "2016APV"]:
        cutAnyFilter = (
            ev.Flag_goodVertices
            and ev.Flag_globalSuperTightHalo2016Filter
            and ev.Flag_HBHENoiseFilter
            and ev.Flag_HBHENoiseIsoFilter
            and ev.Flag_EcalDeadCellTriggerPrimitiveFilter
            and ev.Flag_BadPFMuonFilter
            and ev.Flag_BadPFMuonDzFilter
            and ev.Flag_eeBadScFilter
        )
    else:
        raise ValueError("Unsupported era for quality filters")

    return cutAnyFilter

# Create the actual histograms for saving
histos = {}

if etaOption == "Eta":
    for eta_bin in eta_bins:
        for var in histBins:
            if var == "lep1pt":
                # Use the custom bin edges for lep1pt
                histos[f"{eta_bin}_{var}_num"] = ROOT.TH1F(f"{eta_bin}_{var}_num", f"{eta_bin}_{var}_num", len(histBins[var]) - 1, np.array(histBins[var], dtype=np.float32))
                histos[f"{eta_bin}_{var}_den"] = ROOT.TH1F(f"{eta_bin}_{var}_den", f"{eta_bin}_{var}_den", len(histBins[var]) - 1, np.array(histBins[var], dtype=np.float32))
            else:
                histos[f"{eta_bin}_{var}_num"] = ROOT.TH1F(f"{eta_bin}_{var}_num", f"{eta_bin}_{var}_num", histBins[var][0], histBins[var][1], histBins[var][2])
                histos[f"{eta_bin}_{var}_den"] = ROOT.TH1F(f"{eta_bin}_{var}_den", f"{eta_bin}_{var}_den", histBins[var][0], histBins[var][1], histBins[var][2])
else:
    for var in histBins:
        x_label = pt_label if var == "lep1pt" else var
        if var == "lep1pt":
            # Use the custom bin edges for lep1pt
            histos[var + "_num"] = ROOT.TH1F(var + "_num", var + "_num", len(histBins[var]) - 1, lep1pt_bin_edges)
            histos[var + "_den"] = ROOT.TH1F(var + "_den", var + "_den", len(histBins[var]) - 1, lep1pt_bin_edges)
        else:
            histos[var + "_num"] = ROOT.TH1F(var + "_num", var + "_num", histBins[var][0], histBins[var][1], histBins[var][2])
            histos[var + "_den"] = ROOT.TH1F(var + "_den", var + "_den", histBins[var][0], histBins[var][1], histBins[var][2])

def passes_lepton_cuts(ev, lepton, leptonIndex):
    if lepton == "Muon":
        return (
            ev.Muon_tightId[leptonIndex]
            and abs(ev.Muon_eta[leptonIndex]) < 2.4
            and abs(ev.Muon_dz[leptonIndex]) <= 0.05
            and abs(ev.Muon_dxy[leptonIndex]) <= 0.02
            and ev.Muon_pfIsoId[leptonIndex] >= 5)
    else:
        return (
            ev.Electron_cutBased[leptonIndex] >= 2
            and ev.Electron_mvaFall17V2Iso_WP80[leptonIndex]
            and abs(ev.Electron_dxy[leptonIndex]) < 0.05 + 0.05 * (abs(ev.Electron_eta[leptonIndex]) > 1.479)
            and abs(ev.Electron_dz[leptonIndex]) < 0.10 + 0.10 * (abs(ev.Electron_eta[leptonIndex]) > 1.479)
            and ((abs(ev.Electron_eta[leptonIndex]) < 1.444) or (abs(ev.Electron_eta[leptonIndex]) > 1.566))
            and abs(ev.Electron_eta[leptonIndex]) < 2.5)

def passes_jet_cuts(ev, jetIndex):
        return (
            ev.Jet_pt[jetIndex] > 60)

#calculates dphi between leading lepton and leading jet
def dPhi(obj_phi, jet_phi):
    dphi = obj_phi - jet_phi
    dphi = np.arccos(np.cos(dphi))
    return dphi

def deltaR(eta1, phi1, eta2, phi2):
    """Compute the deltaR between two particles."""
    deta = eta1 - eta2
    dphi = phi1 - phi2
    dphi = np.arctan2(np.sin(dphi), np.cos(dphi))  # Wrap angle between -pi and pi
    return np.sqrt(deta**2 + dphi**2)

def find_leading_lepton(leptons_pt, leptons_eta, leptons_phi):
    """Find the leading lepton (highest pT)."""
    max_pt_index = np.argmax(leptons_pt)
    return leptons_pt[max_pt_index], leptons_eta[max_pt_index], leptons_phi[max_pt_index]

def is_leading_lepton_matched(trigObj_eta, trigObj_phi, trigObj_id, trigObj_filterBits, lepton_eta, lepton_phi, lepton_type):
    """Match the leading lepton to trigger objects."""
    matched = False
    for i in range(len(trigObj_id)):
        if lepton_type == "Electron" and abs(trigObj_id[i]) == 11:
            if ((trigObj_filterBits[i] & 2) == 2) or ((trigObj_filterBits[i] & 2048) == 2048) or ((trigObj_filterBits[i] & 8192) == 8192):
                dR = deltaR(lepton_eta, lepton_phi, trigObj_eta[i], trigObj_phi[i])
                if dR < 0.1:
                    matched = True
                    break
        elif lepton_type == "Muon" and abs(trigObj_id[i]) == 13:
            if (((trigObj_filterBits[i] & 2) == 2) and ((trigObj_filterBits[i] & 8) == 8)) or ((trigObj_filterBits[i] & 1024) == 1024):
                dR = deltaR(lepton_eta, lepton_phi, trigObj_eta[i], trigObj_phi[i])
                if dR < 0.1:
                    matched = True
                    break
    return matched

def get_eta_bin(eta):
    for i, (low, high) in enumerate(eta_ranges):
        if low <= abs(eta) < high:
            return eta_bins[i]
    return None

# Now loop over the events
print("Starting %s" % inputFiles)

inF = 0
nF  = len (inputFiles)

for iFile in inputFiles:
     inF += 1
     print("Starting file %i/%i, %s" % (inF, nF, iFile))
     tf = ROOT.TFile(iFile, "READ")
     events = tf.Get("Events")

     # Event counter
     iEv = 0
     # Total number of events
     nEv = events.GetEntries()

     muon_below27 = 0
     electron_below32 = 0

     # Apply reference cuts for data
     if not passRefCut(events, era):
         continue
         
     for ev in events:
         if iEv % 1000 == 0:
             print("%i/%i events in file done" % (iEv, nEv))
         iEv += 1
#         if(iEv % 1000 == 0): break


         passHLT = False
         # Check if we pass numerator
         for hltpath in hlt:
             if getattr(ev, hltpath, False): passHLT = True

         # Find the lepton with the highest pT
         highest_pt = -1
         highest_pt_lepton_index = -1

         if not getattr(ev, "n" + lepton) > 0:
             continue

         leptons_pt = np.array(getattr(ev, lepton + "_pt"))
         leptons_eta = np.array(getattr(ev, lepton + "_eta"))
         leptons_phi = np.array(getattr(ev, lepton + "_phi"))

         for leptonIndex in range(getattr(ev, "n" + lepton)):
             if passes_lepton_cuts(ev, lepton, leptonIndex):
                 if getattr(ev, lepton + "_pt")[leptonIndex] > highest_pt:
                     highest_pt = getattr(ev, lepton + "_pt")[leptonIndex]
                     highest_pt_lepton_index = leptonIndex

         # If no valid lepton found, continue to the next event
         if highest_pt_lepton_index == -1:
             continue

         #finds leading jet
         highest_jet_pt = -1
         highest_jet_pt_index = -1
         for jetIndex in range(ev.nJet):
             if passes_jet_cuts(ev, jetIndex):
                 if ev.Jet_pt[jetIndex] > highest_jet_pt:
                     highest_jet_pt = ev.Jet_pt[jetIndex]
                     highest_jet_pt_index = jetIndex

         if highest_jet_pt_index == -1 or highest_jet_pt_index == -1:
            continue

         #calculate dphi between jet and lepton and veto it <1.5
         lepton_phi = getattr(ev, lepton + "_phi")[highest_pt_lepton_index]
         jet_phi = ev.Jet_phi[highest_jet_pt_index]
         jet_eta = ev.Jet_eta[highest_jet_pt_index]
         dphi_lepJet = dPhi(lepton_phi, jet_phi)
         deltaR_lepjet = deltaR(leptons_eta, leptons_phi, jet_eta, jet_phi)
         if np.any(deltaR_lepjet < 0.5):
             continue
         if dphi_lepJet < 1.5:
            continue

         #MET cuts
         passmetCut = ev.MET_pt >= offlineCuts["MET"]
         jetIndex = highest_jet_pt_index
         leptonIndex = highest_pt_lepton_index

         lepton_eta = getattr(ev, lepton + "_eta")[highest_pt_lepton_index]
         lepton_phi = getattr(ev, lepton + "_phi")[highest_pt_lepton_index]

         # Trigger object data
         trigObj_eta = np.array(ev.TrigObj_eta)
         trigObj_phi = np.array(ev.TrigObj_phi)
         trigObj_id = np.array(ev.TrigObj_id)
         trigObj_filterBits = np.array(ev.TrigObj_filterBits)

         # Match leading lepton to trigger objects
         lepton_matched = is_leading_lepton_matched(trigObj_eta, trigObj_phi, trigObj_id, trigObj_filterBits, lepton_eta, lepton_phi, lepton)

         if lepton == "Muon" and highest_pt < 27:
             muon_below27 += 1
         elif lepton == "Electron" and highest_pt < 32:
             electron_below32 += 1


         # Variables you want to study
         passlepCut = getattr(ev, lepton + "_pt")[leptonIndex] >= offlineCuts["lep1pt"]
         dphi = ((getattr(ev, lepton + "_phi")[leptonIndex]) - ev.MET_phi)
         mT = (2 * (getattr(ev, lepton + "_pt")[leptonIndex]) *  (ev.MET_pt) * (1 - np.cos(dphi))) ** 0.5
         passmtCut = 30 < mT  < 130

         # Then save denominator and numerator
         if etaOption == "Eta":
            eta_bin = get_eta_bin(lepton_eta)
            if eta_bin is None:
                continue

            for var in histBins:
                passDen = False
                fillvar = None
                if var == "lep1pt":
                    passDen = passes_lepton_cuts(ev, lepton, leptonIndex) and passmetCut and passmtCut
                    fillvar = getattr(ev, lepton + "_pt")[leptonIndex]
                elif var == "MET":
                    passDen = passes_lepton_cuts(ev, lepton, leptonIndex) and passlepCut and passmtCut
                    fillvar = ev.MET_pt
                elif var == "mT":
                    passDen = passes_lepton_cuts(ev, lepton, leptonIndex) and passlepCut and passmetCut
                    fillvar = mT
                elif var == "lep1phi":
                    passDen = passes_lepton_cuts(ev, lepton, leptonIndex) and passmetCut and passmtCut and passlepCut
                    fillvar = getattr(ev, lepton + "_phi")[leptonIndex]
                elif var == "MET phi":
                    passDen = passes_lepton_cuts(ev, lepton, leptonIndex) and passmetCut and passmtCut and passlepCut
                    fillvar = ev.MET_phi
                if passDen and fillvar is not None:
                    histos[f"{eta_bin}_{var}_den"].Fill(fillvar)
                    if passHLT and lepton_matched and refhlt:
                        histos[f"{eta_bin}_{var}_num"].Fill(fillvar)
         else:
            for var in histBins:
                passDen = False
                fillvar = None
                if var == "lep1pt":
                    passDen = passes_lepton_cuts(ev, lepton, leptonIndex) and passmetCut and passmtCut
                    fillvar = getattr(ev, lepton + "_pt")[leptonIndex]
                elif var == "MET":
                    passDen = passes_lepton_cuts(ev, lepton, leptonIndex) and passlepCut and passmtCut
                    fillvar = ev.MET_pt
                elif var == "mT":
                    passDen = passes_lepton_cuts(ev, lepton, leptonIndex) and passlepCut and passmetCut
                    fillvar = mT
                elif var == "lep1phi":
                    passDen = passes_lepton_cuts(ev, lepton, leptonIndex) and passmetCut and passmtCut and passlepCut
                    fillvar = getattr(ev, lepton + "_phi")[leptonIndex]
                elif var == "MET phi":
                    passDen = passes_lepton_cuts(ev, lepton, leptonIndex) and passmetCut and passmtCut and passlepCut
                    fillvar = ev.MET_phi
                if passDen and fillvar is not None:
                    histos[var + "_den"].Fill(fillvar)
                    if passHLT and lepton_matched and refhlt:
                        histos[var + "_num"].Fill(fillvar)

     tf.Close()

for var in histBins:
    if etaOption == "Eta":
        for eta_bin in eta_bins:
            print(f"Number of events in the denominator for {eta_bin}_{var}: {histos[f'{eta_bin}_{var}_den'].GetEntries()}")
            print(f"Number of events passing in the numerator for {eta_bin}_{var}: {histos[f'{eta_bin}_{var}_num'].GetEntries()}")
    else:
        print(f"Number of events in the denominator for {var}: {histos[var + '_den'].GetEntries()}")
        print(f"Number of events passing in the numerator for {var}: {histos[var + '_num'].GetEntries()}")

outF = ROOT.TFile(outputHistos, "RECREATE")
for h in histos:
    histos[h].Write()
outF.Close()

print("Processing complete. Output written to %s" % outputHistos)

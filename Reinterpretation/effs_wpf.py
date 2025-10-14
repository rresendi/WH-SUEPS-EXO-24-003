import ROOT
import os
import matplotlib.pyplot as plt
import numpy as np
import mplhep as hep
from array import array
import sys
import csv
from math import isnan

# Debugging toggle
DEBUG = False

def dbg(*args, **kwargs):
    if DEBUG:
        print(*args, **kwargs)

# Confidence level and alpha for Clopper-Pearson
CL = 0.683
alpha = 1.0 - CL  # TEfficiency.ClopperPearson expects alpha = 1 - CL

# DeltaR helper function
def deltaR(eta1, phi1, eta2, phi2):
    dphi = np.abs(phi1 - phi2)
    if dphi > np.pi:
        dphi = 2*np.pi - dphi
    deta = eta1 - eta2
    return np.sqrt(deta**2 + dphi**2)

# Helper to check bin membership
def in_bin(pt, aeta, pt_lo, pt_hi, eta_lo, eta_hi):
    return (pt_lo <= pt < pt_hi) and (eta_lo <= aeta < eta_hi)

def trace_mother_chain(gen_idx, mother_indices, pdg_ids, target_pdg=None):
    """Trace the GenPart mother chain and flag if a target PDG appears."""
    chain_indices = []
    chain_pdg_ids = []
    found_target = target_pdg is None
    current = gen_idx
    visited = set()

    while True:
        if current < 0 or current >= len(mother_indices):
            break
        mother_idx = int(mother_indices[current])
        if mother_idx < 0 or mother_idx >= len(pdg_ids):
            break
        if mother_idx in visited:
            break
        visited.add(mother_idx)
        pdg = int(pdg_ids[mother_idx])
        chain_indices.append(mother_idx)
        chain_pdg_ids.append(pdg)
        if target_pdg is not None and abs(pdg) == target_pdg:
            found_target = True
        current = mother_idx

    return found_target, chain_indices, chain_pdg_ids

def main():

    input_dir = sys.argv[1] # Path to input files
    samp = sys.argv[2] # dy, jpsi, or zprime
    obj = "muon"

    #Gather input files and initialize file and event indices
    inputFiles = [os.path.join(input_dir, f) for f in os.listdir(input_dir) if f.endswith(".root")]
    nF = len(inputFiles)

    # Set up 2D histograms with binning used in current nominal Delphes card
    eta_edges = array('d', [0, 0.9, 1.2, 2.1, 2.4])

    if samp == "jpsi":
        pt_edges = array('d', [5.0, 10.0, 15.0])
    elif samp == "dy":
        pt_edges = array('d', [15, 100, 1000])
    elif samp == "zprime":
        pt_edges = array('d', [1000, 3000])

    if samp == "jpsi":
        resonance_mass = 3.0969     # J/psi
    elif samp == "dy":
        resonance_mass = 91.1876    # Z
    elif samp == "zprime":
        resonance_mass = 2500.0     # z'
    else:
        resonance_mass = 0.0

    target_lookup = {
        "jpsi": 443,
        "dy": 23,
        "zprime": 32,
    }
    target_label_lookup = {
        "jpsi": "J/psi",
        "dy": "Z",
        "zprime": "Z'",
    }
    samp_key = samp.lower()
    target_pdg = target_lookup.get(samp_key)
    target_label = target_label_lookup.get(samp_key, "target resonance")

    # Global assigned and unassigned gen leptons
    assigned_gen, unassigned_gen = [], []

    # Event-level rows for plotting (reco-level OSSF pair): mass and leading/subleading kinematics
    event_rows = []  # each: {mass, lead_pt, lead_eta, lead_phi, sub_pt, sub_eta, sub_phi}

    # Global assigned and unassigned pf candidates
    assigned_pf, unassigned_pf = [], []

    # Store momentum resolution for gen→PF matched objects (one row per match)
    # Each entry: {gen_eta, pf_eta, gen_pt, pf_pt, resolution}
    mom_resolutions = []

    # Track GenPart ancestry information
    gen_mother_records = []
    record_lookup = {}
    pf_mismatch_records = []
    total_gen_selected = 0
    total_gen_from_target = 0

    def append_gen_mother_records(gen_collection, best_indices, file_index, file_name, event_index):
        for obj in gen_collection:
            chain_idx = [int(idx) for idx in obj.get("mother_chain_indices", [])]
            chain_pdg = [int(pid) for pid in obj.get("mother_chain_pdgIds", [])]
            mother_pdg = chain_pdg[0] if chain_pdg else None
            mother_non_target = None
            if mother_pdg is not None and (target_pdg is None or abs(mother_pdg) != target_pdg):
                mother_non_target = mother_pdg
            key = (file_index, event_index, int(obj["index"]))
            gen_mother_records.append({
                "file_index": file_index,
                "file_name": file_name,
                "event_index": event_index,
                "gen_index": int(obj["index"]),
                "pdgId": int(obj["pdgId"]),
                "pt": float(obj["pt"]),
                "eta": float(obj["eta"]),
                "phi": float(obj["phi"]),
                "mass": float(obj["mass"]),
                "from_target": bool(obj.get("from_target", False)),
                "in_best_pair": obj["index"] in best_indices,
                "mother_chain_indices": chain_idx,
                "mother_chain_pdgIds": chain_pdg,
                "mother_pdgId_non_target": mother_non_target,
                "pf_mismatch_pdgId": None,
            })
            record_lookup[key] = len(gen_mother_records) - 1

    # Start looping through files in the directory
    for z, iFile in enumerate(inputFiles):
        print(f"Processing file {z}/{nF}: {iFile}")
        file_basename = os.path.basename(iFile)

        # Open ROOT file and get Events tree safely
        try:
            tf = ROOT.TFile.Open(iFile, "READ")
            if not tf or tf.IsZombie():
                raise OSError(f"File {iFile} is unreadable or corrupted.")
            events = tf.Get("Events")
            if not events:
                raise OSError(f"No 'Events' tree found in {iFile}.")
        except OSError as e:
            print(f"Error opening file {iFile}: {e} Skipping.")
            if 'tf' in locals() and tf:  # Close file if it was opened
                tf.Close()
            continue

        # Start looping through the events in a given file
        for iEv, ev in enumerate(events):
            #For debugging
            if(iEv == 10000):
                break

            # For event header debug
            dbg(f"[Event {iEv}] processing...")

            # Extract relevant branches
            nGenPart = ev.nGenPart
            pdgIds = ev.GenPart_pdgId
            pts = ev.GenPart_pt
            etas = ev.GenPart_eta
            phis = ev.GenPart_phi
            status = ev.GenPart_status
            masses_branch = ev.GenPart_mass
            mother_indices = ev.GenPart_genPartIdxMother

            # Create gen particle collection with pT-range requirement from pt_edges
            gen_objs = []
            pt_min, pt_max = pt_edges[0], pt_edges[-1]
            if obj == "muon":
                gen_pdgid, gen_eta = 13, 2.4
            else:
                gen_pdgid, gen_eta = 11, 2.5

            for i in range(nGenPart):
                # Apply object-level selections, including pT range
                if (
                    abs(pdgIds[i]) == gen_pdgid
                    and abs(etas[i]) < gen_eta
                    and status[i] == 1
                    and (pts[i] >= pt_min and pts[i] < pt_max)
                ):
                    from_target = False
                    chain_indices, chain_pdgids = [], []
                    mother_idx = int(mother_indices[i]) if i < len(mother_indices) else -1
                    if 0 <= mother_idx < len(pdgIds):
                        mother_pdg = int(pdgIds[mother_idx])
                        if abs(mother_pdg) == target_pdg:
                            from_target = True
                    # optional fallback
                    if not from_target:
                        from_target, chain_indices, chain_pdgids = trace_mother_chain(i, mother_indices, pdgIds, target_pdg)
                    else:
                        from_target, chain_indices, chain_pdgids = True, [mother_idx], [mother_pdg]
                    if not from_target:
                        continue
                    total_gen_selected += 1
                    if from_target:
                        total_gen_from_target += 1
                    gen_objs.append({
                        "pt": pts[i],
                        "eta": etas[i],
                        "phi": phis[i],
                        "mass": masses_branch[i],
                        "pdgId": pdgIds[i],
                        "index": i,
                        "from_target": from_target,
                        "mother_chain_indices": chain_indices,
                        "mother_chain_pdgIds": chain_pdgids
                    })
                        
            # Sort by pt
            gen_objs = sorted(gen_objs, key=lambda x: x["pt"], reverse=True)

            # At least two candidates per event
            if len(gen_objs) < 2:
                append_gen_mother_records(gen_objs, set(), z, file_basename, iEv)
                continue

            best_pair = None
            best_delta = float("inf")

            # Loop through all gen particles for the given flavor
            for i, obj1 in enumerate(gen_objs):
                for obj2 in gen_objs[i+1:]:

                    # Do the OSSF check
                    if obj1["pdgId"] * obj2["pdgId"] >= 0:
                        continue

                    # Calculate resonance invariant mass only using OSSF pairs and apply resonance mass selection
                    p4_1 = ROOT.TLorentzVector()
                    p4_1.SetPtEtaPhiM(obj1["pt"], obj1["eta"], obj1["phi"], obj1["mass"])
                    p4_2 = ROOT.TLorentzVector()
                    p4_2.SetPtEtaPhiM(obj2["pt"], obj2["eta"], obj2["phi"], obj2["mass"])

                    inv_mass = (p4_1 + p4_2).M()
                    # print("gen_inv: ", inv_mass)

                    if samp == "jpsi":
                        if inv_mass < 2.9 or inv_mass > 3.3:
                            continue
                    elif samp == "dy":
                        if inv_mass < 60 or inv_mass > 120:
                            continue
                    elif samp == "zprime":
                        if inv_mass < 2300 or inv_mass > 2700:
                            continue
                        
                    delta = abs(inv_mass - resonance_mass)
                    if delta < best_delta:
                        best_delta = delta
                        best_pair = (obj1, obj2)

            # If no OSSF pair, skip this event
            if best_pair is None:
                append_gen_mother_records(gen_objs, set(), z, file_basename, iEv)
                continue

            best_indices = {best_pair[0]["index"], best_pair[1]["index"]}
            append_gen_mother_records(gen_objs, best_indices, z, file_basename, iEv)

            # Adding particles to the local event and global OSSF gen-muons set
            best_gen_pair_objs = [best_pair[0], best_pair[1]]

            # Build PF candidate list with pion/electron/photon vetoes
            pf_candidates = []
            count = 0
            if obj == "muon":
                veto_pdgids = {211, 11, 22}
            elif obj == "electron":
                veto_pdgids = {211, 13, 22}
            else:
                veto_pdgids = {211, 22}

            for i in range(ev.nPFCands):
                if abs(ev.PFCands_pdgId[i]) in veto_pdgids:
                    continue
                pf_candidates.append({
                    "index": count,
                    "pt": ev.PFCands_trkPt[i],
                    "eta": ev.PFCands_trkEta[i],
                    "phi": ev.PFCands_trkPhi[i],
                    "mass": ev.PFCands_mass[i],
                    "pdgid": ev.PFCands_pdgId[i],
                })
                count += 1

            dbg(f"[Event {iEv}] nPF={len(pf_candidates)}")

            # Build all valid (lepton j, PF i) pairs passing the selection, scored by (dR, dPtRel)
            gen_pairs = []
            for j, gen_obj in enumerate(best_gen_pair_objs):
                gpt = gen_obj["pt"]
                geta = gen_obj["eta"]
                gphi = gen_obj["phi"]
                for i, pf in enumerate(pf_candidates):
                    dR = deltaR(geta, gphi, pf["eta"], pf["phi"])
                    if dR >= 0.05:
                        continue
                    if gpt:
                        dPtRel = abs(gpt - pf["pt"]) / gpt
                    else:
                        dPtRel = float("inf")
                    if obj == "muon" and dPtRel >= 0.3:
                        continue
                    # Score: prioritize smallest dR, then smallest dPtRel
                    gen_pairs.append((dR, dPtRel, j, i))

            # Greedy one-to-one assignment by best (dR, dPtRel)
            gen_pairs.sort()  # sorts by dR first, then dPtRel
            used_pf_indices = set()   # local helper: indices of PF candidates already consumed
            assigned_leptons = {}

            for dR, dPtRel, j, i in gen_pairs:
                if j in assigned_leptons or i in used_pf_indices:
                    continue
                assigned_leptons[j] = {
                    "pf_list_idx": i,
                    "pf_index": pf_candidates[i]["index"],
                    "dR": dR,
                    "dPtRel": dPtRel
                }
                used_pf_indices.add(i)
                if len(assigned_leptons) == 2:
                    break

            # After assignment, record whether each lepton was assigned or unassigned.
            # Append the lepton 4-momentum (pt, eta, phi, mass) to the appropriate global arrays.
            for j, gen_obj in enumerate(best_gen_pair_objs):
                fourvec = {
                    "pt": gen_obj["pt"],
                    "eta": gen_obj["eta"],
                    "phi": gen_obj["phi"],
                    "mass": gen_obj["mass"],
                    "from_target": gen_obj.get("from_target", False)
                }
                if j in assigned_leptons:
                    assigned_gen.append(fourvec)
                else:
                    unassigned_gen.append(fourvec)

            # Store assigned PfCands (keep dicts so we can access by keys)
            pf_matched_to_gen = []
            expected_abs_pdg = 13

            for j, gen_obj in enumerate(best_gen_pair_objs):
                if j in assigned_leptons:
                    gpt = gen_obj["pt"]
                    geta = gen_obj["eta"]
                    gphi = gen_obj["phi"]
                    gmass = gen_obj["mass"]
                    pf_list_idx = assigned_leptons[j]["pf_list_idx"]
                    pf = pf_candidates[pf_list_idx]
                    pf_report_index = assigned_leptons[j]["pf_index"]
                    pf_pdg = int(pf["pdgid"])
                    pf_matched_to_gen.append({
                        "pt": pf["pt"],
                        "eta": pf["eta"],
                        "phi": pf["phi"],
                        "mass": pf["mass"],
                        "pf_index": pf_report_index,
                        # carry corresponding gen information for resolution
                        "gen_pt": gpt,
                        "gen_eta": geta,
                        "gen_phi": gphi,
                        "gen_mass": gmass,
                        "from_target": gen_obj.get("from_target", False),
                        "pf_pdgId": pf_pdg
                    })
                    if abs(pf_pdg) != expected_abs_pdg:
                        pf_mismatch_records.append({
                            "file_index": z,
                            "file_name": file_basename,
                            "event_index": iEv,
                            "gen_index": int(gen_obj["index"]),
                            "gen_pdgId": int(gen_obj["pdgId"]),
                            "pf_index": pf_report_index,
                            "pf_pdgId": pf_pdg,
                            "pf_pt": float(pf["pt"]),
                            "expected_abs_pdgId": expected_abs_pdg
                        })
                        lookup_key = (z, iEv, int(gen_obj["index"]))
                        rec_idx = record_lookup.get(lookup_key)
                        if rec_idx is not None:
                            gen_mother_records[rec_idx]["pf_mismatch_pdgId"] = pf_pdg

            # Initialize reconstructed objects for second efficiency
            reco_objs = []
            for i in range(ev.nMuon):
                reco_objs.append({
                "pt": ev.Muon_pt[i],
                "eta": ev.Muon_eta[i],
                "phi": ev.Muon_phi[i],
                "mass": ev.Muon_mass[i]
            })

            dbg(f"[Event {iEv}] nReco={len(reco_objs)}  nPF_matched_to_gen={len(pf_matched_to_gen)}")

            # Compute momentum resolution for gen→PF matched objects
            # res = (pf_pt - gen_pt) / gen_pt, for all PF matched to gen
            for pf in pf_matched_to_gen:
                gen_pt = pf["gen_pt"]
                if gen_pt and gen_pt > 0:
                    res = (pf["pt"] - gen_pt) / gen_pt
                    mom_resolutions.append({
                        "gen_eta": float(pf["gen_eta"]),
                        "pf_eta": float(pf["eta"]),
                        "gen_pt": float(gen_pt),
                        "pf_pt": float(pf["pt"]),
                        "resolution": float(res)
                    })

            # Build all valid (PF j, reco i) pairs passing the selection, scored by (dR, dPtRel)
            reco_pairs = []
            for j, reco in enumerate(reco_objs):
                for i, pf in enumerate(pf_matched_to_gen):
                    dR = deltaR(reco["eta"], reco["phi"], pf["eta"], pf["phi"])
                    if dR >= 0.05:
                        continue
                    dPtRel = abs(reco["pt"] - pf["pt"]) / reco["pt"]
                    if obj == "muon" and dPtRel >= 0.3:
                        continue
                    reco_pairs.append((dR, dPtRel, j, i))
            dbg(f"[Event {iEv}] reco_pairs (candidate links) = {len(reco_pairs)}")

            # Greedy one-to-one assignment by best (dR, dPtRel)
            reco_pairs.sort()  # sorts by dR first, then dPtRel
            assigned_reco = set()             # set of reco indices already matched # dR match between this and best_gen_pair_objs to check which is closest to reco and use that for
            assigned_pf_map = {}              # map: reco j  -> pf i (index in pf_matched_to_gen)

            for dR, dPtRel, j, i in reco_pairs:
                if j in assigned_pf_map or i in assigned_reco:
                    continue
                assigned_pf_map[j] = i
                assigned_reco.add(i)
                if len(assigned_pf_map) == len(pf_matched_to_gen):
                    break

            dbg(f"[Event {iEv}] assigned reco→PF pairs = {len(assigned_pf_map)} of {len(pf_matched_to_gen)} PF matched-to-gen")

            # After assignment, record whether each PF (matched-to-gen) was assigned to a reco object or not.
            # Append the PF 4-momentum (pt, eta, phi, mass) to the appropriate global arrays.
            assigned_pf_indices = set(assigned_pf_map.values())
            for i, pf in enumerate(pf_matched_to_gen):
                fourvec_pf = (pf["pt"], pf["eta"], pf["phi"], pf["mass"]) 
                if i in assigned_pf_indices:
                    assigned_pf.append(fourvec_pf)
                else:
                    unassigned_pf.append(fourvec_pf)
            assigned_now = sum(1 for i in assigned_pf_map.values())
            unassigned_now = len(pf_matched_to_gen) - assigned_now
            dbg(f"[Event {iEv}] PF→reco: assigned={assigned_now}, unassigned={unassigned_now}")

            # Compute reconstructed invariant mass if two matched reco leptons are available
            matched_reco = []
            for reco_idx, pf_idx in assigned_pf_map.items():
                if 0 <= reco_idx < len(reco_objs) and 0 <= pf_idx < len(pf_matched_to_gen):
                    reco = reco_objs[reco_idx]
                    matched_reco.append({
                        "pt": reco["pt"],
                        "eta": reco["eta"],
                        "phi": reco["phi"],
                        "mass": reco["mass"]
                    })

            if len(matched_reco) >= 2:
                matched_reco.sort(key=lambda x: x["pt"], reverse=True)
                reco_p4_1 = ROOT.TLorentzVector()
                reco_p4_1.SetPtEtaPhiM(
                    matched_reco[0]["pt"], matched_reco[0]["eta"], matched_reco[0]["phi"], matched_reco[0]["mass"]
                )
                reco_p4_2 = ROOT.TLorentzVector()
                reco_p4_2.SetPtEtaPhiM(
                    matched_reco[1]["pt"], matched_reco[1]["eta"], matched_reco[1]["phi"], matched_reco[1]["mass"]
                )
                event_rows.append({
                    "mass": float((reco_p4_1 + reco_p4_2).M()),
                    "lead_pt": float(matched_reco[0]["pt"]),
                    "lead_eta": float(matched_reco[0]["eta"]),
                    "lead_phi": float(matched_reco[0]["phi"]),
                    "sub_pt": float(matched_reco[1]["pt"]),
                    "sub_eta": float(matched_reco[1]["eta"]),
                    "sub_phi": float(matched_reco[1]["phi"]),
                })
        tf.Close()

    print(f"Finished processing {nF} files.")
    print(f"Total selected gen {obj}s: {total_gen_selected}")
    if target_pdg is not None:
        print(
            f"Total selected gen {obj}s with {target_label} ancestor (|PDG|={target_pdg}): {total_gen_from_target}"
        )
    else:
        print(f"No mother PDG requirement applied; denominator uses all selected gen {obj}s.")

    # Define bin lists
    pt_bins = [(pt_edges[i], pt_edges[i+1]) for i in range(len(pt_edges)-1)]
    eta_bins = [(eta_edges[i], eta_edges[i+1]) for i in range(len(eta_edges)-1)]
    
    # Prepare gen containers
    n_eta = len(eta_bins)
    n_pt = len(pt_bins)
    eff_matrix_gen = [[float('nan') for _ in range(n_pt)] for _ in range(n_eta)]
    err_lo_matrix_gen = [[0.0 for _ in range(n_pt)] for _ in range(n_eta)]
    err_up_matrix_gen = [[0.0 for _ in range(n_pt)] for _ in range(n_eta)]
    assigned_counts_gen = [[0 for _ in range(n_pt)] for _ in range(n_eta)]
    total_counts_gen = [[0 for _ in range(n_pt)] for _ in range(n_eta)]

    # Prepare pf containers
    eff_matrix_pf = [[float('nan') for _ in range(n_pt)] for _ in range(n_eta)]
    err_lo_matrix_pf = [[0.0 for _ in range(n_pt)] for _ in range(n_eta)]
    err_up_matrix_pf = [[0.0 for _ in range(n_pt)] for _ in range(n_eta)]
    assigned_counts_pf = [[0 for _ in range(n_pt)] for _ in range(n_eta)]
    total_counts_pf = [[0 for _ in range(n_pt)] for _ in range(n_eta)]

    # Debug: inspect PF bin denominators/totals
    dbg(f"[Global] total assigned PF 4vecs = {len(assigned_pf)}")
    dbg(f"[Global] total unassigned PF 4vecs = {len(unassigned_pf)}")
    
    # Count assigned gen
    for entry in assigned_gen:
        if not entry.get("from_target", False):
            continue
        pt = entry["pt"]
        aeta = abs(entry["eta"])
        for ie, (eta_lo, eta_hi) in enumerate(eta_bins):
            if not (eta_lo <= aeta < eta_hi):
                continue
            for ip, (pt_lo, pt_hi) in enumerate(pt_bins):
                if pt_lo <= pt < pt_hi:
                    assigned_counts_gen[ie][ip] += 1
                    total_counts_gen[ie][ip] += 1
                    break
    
    # Count unassigned gen
    for entry in unassigned_gen:
        if not entry.get("from_target", False):
            continue
        pt = entry["pt"]
        aeta = abs(entry["eta"])
        for ie, (eta_lo, eta_hi) in enumerate(eta_bins):
            if not (eta_lo <= aeta < eta_hi):
                continue
            for ip, (pt_lo, pt_hi) in enumerate(pt_bins):
                if pt_lo <= pt < pt_hi:
                    total_counts_gen[ie][ip] += 1
                    break

    # Count assigned pf
    for (pt, eta, _, _) in assigned_pf:
        aeta = abs(eta)
        for ie, (eta_lo, eta_hi) in enumerate(eta_bins):
            if not (eta_lo <= aeta < eta_hi):
                continue
            for ip, (pt_lo, pt_hi) in enumerate(pt_bins):
                if pt_lo <= pt < pt_hi:
                    assigned_counts_pf[ie][ip] += 1
                    total_counts_pf[ie][ip] += 1
                    break

    # Count unassigned pf
    for (pt, eta, _, _) in unassigned_pf:
        aeta = abs(eta)
        for ie, (eta_lo, eta_hi) in enumerate(eta_bins):
            if not (eta_lo <= aeta < eta_hi):
                continue
            for ip, (pt_lo, pt_hi) in enumerate(pt_bins):
                if pt_lo <= pt < pt_hi:
                    total_counts_pf[ie][ip] += 1
                    break

    # Debug: list PF bins with zero denominator
    for ie, (eta_lo, eta_hi) in enumerate(eta_bins):
        for ip, (pt_lo, pt_hi) in enumerate(pt_bins):
            if total_counts_pf[ie][ip] == 0:
                dbg(f"[PF Bin Empty] eta=[{eta_lo},{eta_hi}) pt=[{pt_lo},{pt_hi}) assigned={assigned_counts_pf[ie][ip]} total={total_counts_pf[ie][ip]}")
    
    # Compute gen efficiency and Clopper-Pearson bounds per bin (using ROOT.TEfficiency static method)
    for ie in range(n_eta):
        for ip in range(n_pt):
            k = assigned_counts_gen[ie][ip]
            n = total_counts_gen[ie][ip]
            if n > 0:
                eff = k / n
                low = ROOT.TEfficiency.ClopperPearson(n, k, alpha, False)  # lower bound
                up  = ROOT.TEfficiency.ClopperPearson(n, k, alpha, True)   # upper bound
                eff_matrix_gen[ie][ip] = eff
                err_lo_matrix_gen[ie][ip] = eff - low
                err_up_matrix_gen[ie][ip] = up - eff
            else:
                eff_matrix_gen[ie][ip] = float('nan')
                err_lo_matrix_gen[ie][ip] = 0.0
                err_up_matrix_gen[ie][ip] = 0.0

    # Compute pf efficiency
    for ie in range(n_eta):
        for ip in range(n_pt):
            k = assigned_counts_pf[ie][ip]
            n = total_counts_pf[ie][ip]
            if n > 0:
                eff = k / n
                low = ROOT.TEfficiency.ClopperPearson(n, k, alpha, False)  # lower bound
                up  = ROOT.TEfficiency.ClopperPearson(n, k, alpha, True)   # upper bound
                eff_matrix_pf[ie][ip] = eff
                err_lo_matrix_pf[ie][ip] = eff - low
                err_up_matrix_pf[ie][ip] = up - eff
            else:
                dbg(f"[PF Eff NaN] eta bin {ie}, pt bin {ip} has total=0; leaving eff as NaN")
                eff_matrix_pf[ie][ip] = float('nan')
                err_lo_matrix_pf[ie][ip] = 0.0
                err_up_matrix_pf[ie][ip] = 0.0

    # Human-readable labels
    pt_labels = [f"[{lo:.3g},{hi:.3g})" for (lo, hi) in pt_bins]
    eta_labels = [f"[{lo:.3g},{hi:.3g})" for (lo, hi) in eta_bins]
    
    # Provide a tidy (long-form) CSV with counts and errors per bin
    with open(f"{samp}_{obj}_m1_efficiency.csv", "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["eta_lo", "eta_hi", "pt_lo", "pt_hi", "gen_assigned", "gen_total", "gen_eff", "gen_err_low", "gen_err_up", "pf_assigned", "pf_total", "pf_eff", "pf_err_low", "pf_err_up"])
        for ie, (eta_lo, eta_hi) in enumerate(eta_bins):
            for ip, (pt_lo, pt_hi) in enumerate(pt_bins):
                # Gen eff info
                k = assigned_counts_gen[ie][ip]
                n = total_counts_gen[ie][ip]
                gen_eff = eff_matrix_gen[ie][ip]
                err_lo_gen = err_lo_matrix_gen[ie][ip]
                err_up_gen = err_up_matrix_gen[ie][ip]
                # Pf eff info
                l = assigned_counts_pf[ie][ip]
                m = total_counts_pf[ie][ip]
                pf_eff = eff_matrix_pf[ie][ip]
                err_lo_pf = err_lo_matrix_pf[ie][ip]
                err_up_pf = err_up_matrix_pf[ie][ip]
                # Write NaN as empty string for clarity
                gen_eff_str = "" if isnan(gen_eff) else f"{gen_eff:.6f}"
                pf_eff_str = "" if isnan(pf_eff) else f"{pf_eff:.6f}"
                writer.writerow([
                    eta_lo, eta_hi, pt_lo, pt_hi, k, n, gen_eff_str, f"{err_lo_gen:.6f}", f"{err_up_gen:.6f}", l, m, pf_eff_str, f"{err_lo_pf:.6f}",  f"{err_up_pf:.6f}"
                ])
    
    print(f"Saved assignment efficiency CSVs: {samp}_{obj}_m1_efficiency.csv")

    with open(f"{samp}_{obj}_m1_gen_ancestry.csv", "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow([
            "file_index",
            "file_name",
            "event_index",
            "gen_index",
            "pdgId",
            "pt",
            "eta",
            "phi",
            "mass",
            "from_target",
            "in_best_pair",
            "mother_chain_indices",
            "mother_chain_pdgIds",
            "mother_pdgId_if_not_target",
            "pf_mismatch_pdgId",
        ])
        for rec in gen_mother_records:
            chain_idx_str = "->".join(str(idx) for idx in rec["mother_chain_indices"])
            chain_pdg_str = "->".join(str(pid) for pid in rec["mother_chain_pdgIds"])
            mother_non_target = rec.get("mother_pdgId_non_target")
            mother_non_target_str = "" if mother_non_target is None else str(mother_non_target)
            pf_mismatch_val = rec.get("pf_mismatch_pdgId")
            pf_mismatch_str = "" if pf_mismatch_val is None else str(pf_mismatch_val)
            writer.writerow([
                rec["file_index"],
                rec["file_name"],
                rec["event_index"],
                rec["gen_index"],
                rec["pdgId"],
                f"{rec['pt']:.6f}",
                f"{rec['eta']:.6f}",
                f"{rec['phi']:.6f}",
                f"{rec['mass']:.6f}",
                1 if rec["from_target"] else 0,
                1 if rec["in_best_pair"] else 0,
                chain_idx_str,
                chain_pdg_str,
                mother_non_target_str,
                pf_mismatch_str,
            ])
    print(f"Saved gen ancestry CSV: {samp}_{obj}_m1_gen_ancestry.csv")

    mismatch_filename = f"{samp}_{obj}_m1_pf_mismatches.csv"
    if pf_mismatch_records:
        with open(mismatch_filename, "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow([
                "file_index",
                "file_name",
                "event_index",
                "gen_index",
                "gen_pdgId",
                "pf_index",
                "pf_pdgId",
                "pf_pt",
                "expected_abs_pdgId",
            ])
            for rec in pf_mismatch_records:
                writer.writerow([
                    rec["file_index"],
                    rec["file_name"],
                    rec["event_index"],
                    rec["gen_index"],
                    rec["gen_pdgId"],
                    rec["pf_index"],
                    rec["pf_pdgId"],
                    f"{rec['pf_pt']:.6f}",
                    rec["expected_abs_pdgId"],
                ])
        print(f"Saved PF mismatch CSV: {mismatch_filename} ({len(pf_mismatch_records)} entries)")
    else:
        print("No PF PDG ID mismatches found for matched leptons.")

    # Write per-event gen-pair CSV for plotting distributions
    with open(f"{samp}_{obj}_m1_kins.csv", "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["mass", "lead_pt", "lead_eta", "lead_phi", "sub_pt", "sub_eta", "sub_phi"])
        for r in event_rows:
            mass_val = r['mass']
            mass_str = "" if isnan(mass_val) else f"{mass_val:.6f}"
            writer.writerow([
                mass_str,
                f"{r['lead_pt']:.6f}", f"{r['lead_eta']:.6f}", f"{r['lead_phi']:.6f}",
                f"{r['sub_pt']:.6f}",  f"{r['sub_eta']:.6f}",  f"{r['sub_phi']:.6f}",
            ])
    print(f"Saved event-level CSV: {samp}_{obj}_m1_kins.csv")

    # Write per-match momentum resolution to its own CSV
    with open(f"{samp}_{obj}_m1_momentum_resolution.csv", "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["gen_eta", "pf_eta", "gen_pt", "pf_pt", "resolution"])
        for entry in mom_resolutions:
            writer.writerow([
                f"{entry['gen_eta']:.6f}",
                f"{entry['pf_eta']:.6f}",
                f"{entry['gen_pt']:.6f}",
                f"{entry['pf_pt']:.6f}",
                f"{entry['resolution']:.6f}"
            ])
    print(f"Saved momentum resolution CSV: {samp}_{obj}_m1_momentum_resolution.csv")

if __name__ == "__main__":
    main()

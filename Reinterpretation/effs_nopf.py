import ROOT
import os
import sys
import csv
import numpy as np

# Confidence level and alpha for Clopper-Pearson
CL = 0.683
alpha = 1.0 - CL  # TEfficiency.ClopperPearson expects alpha = 1 - CL

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
    if len(sys.argv) < 3:
        print("Usage: python muon_tracking_efficiency.py <input_dir> <sample>")
        sys.exit(1)

    input_dir = sys.argv[1]  # Path to input files
    samp = sys.argv[2].lower()  # dy, jpsi, or zprime
    obj = "muon"

    # Gather input files
    inputFiles = [os.path.join(input_dir, f) for f in os.listdir(input_dir) if f.endswith(".root")]
    nF = len(inputFiles)

    # Binning
    eta_edges = [0, 0.9, 1.2, 2.1, 2.4]
    if samp == "jpsi":
        pt_edges = [5.0, 10.0, 15.0]
        resonance_mass = 3.0969
    elif samp == "dy":
        pt_edges = [15, 100, 1000]
        resonance_mass = 91.1876
    elif samp == "zprime":
        pt_edges = [1000, 3000]
        resonance_mass = 2500.0
    else:
        raise ValueError("sample must be one of: jpsi, dy, zprime")

    target_lookup = {"jpsi": 443, "dy": 23, "zprime": 32}
    target_label_lookup = {"jpsi": "J/psi", "dy": "Z", "zprime": "Z'"}
    target_pdg = target_lookup[samp]
    target_label = target_label_lookup[samp]

    # Collections
    event_rows = []  # reco OSSF pair kinematics
    gen_status_records = []  # per-gen status used to compute efficiencies
    reco_momentum_resolutions = []
    gen_mother_records = []
    record_lookup = {}
    reco_mismatch_records = []
    total_gen_selected = 0
    total_gen_from_target = 0

    def append_gen_mother_records(gen_collection, best_indices, file_index, file_name, event_index):
        for obj in gen_collection:
            chain_idx = obj.get("mother_chain_indices", [])
            chain_pdg = obj.get("mother_chain_pdgIds", [])
            mother_pdg = chain_pdg[0] if chain_pdg else None
            mother_non_target = None
            if mother_pdg is not None and abs(mother_pdg) != target_pdg:
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
                "reco_mismatch_pdgId": None,
            })
            record_lookup[key] = len(gen_mother_records) - 1

    # Loop over files
    for z, iFile in enumerate(inputFiles):
        print(f"Processing file {z}/{nF}: {iFile}")
        try:
            tf = ROOT.TFile.Open(iFile, "READ")
            if not tf or tf.IsZombie():
                raise OSError(f"File {iFile} is unreadable or corrupted.")
            events = tf.Get("Events")
            if not events:
                raise OSError(f"No 'Events' tree found in {iFile}.")
        except OSError as e:
            print(f"Error opening file {iFile}: {e} Skipping.")
            if 'tf' in locals() and tf:
                tf.Close()
            continue

        # Events
        for iEv, ev in enumerate(events):
            # Extract relevant branches
            nGenPart = ev.nGenPart
            pdgIds = ev.GenPart_pdgId
            pts = ev.GenPart_pt
            etas = ev.GenPart_eta
            phis = ev.GenPart_phi
            status = ev.GenPart_status
            masses_branch = ev.GenPart_mass
            mother_indices = ev.GenPart_genPartIdxMother

            # Select gen muons from immediate target decays within kinematic range
            gen_objs = []
            pt_min, pt_max = pt_edges[0], pt_edges[-1]
            for i in range(nGenPart):
                if (
                    abs(pdgIds[i]) == 13
                    and abs(etas[i]) < 2.4
                    and status[i] == 1
                    and (pt_min <= pts[i] < pt_max)
                ):
                    mother_idx = int(mother_indices[i]) if i < len(mother_indices) else -1
                    if 0 <= mother_idx < len(pdgIds) and abs(int(pdgIds[mother_idx])) == target_pdg:
                        _, chain_indices, chain_pdgids = trace_mother_chain(i, mother_indices, pdgIds, target_pdg)
                        total_gen_selected += 1
                        total_gen_from_target += 1
                        gen_objs.append({
                            "pt": float(pts[i]),
                            "eta": float(etas[i]),
                            "phi": float(phis[i]),
                            "mass": float(masses_branch[i]),
                            "pdgId": int(pdgIds[i]),
                            "index": i,
                            "from_target": True,
                            "mother_index": mother_idx,
                            "mother_pdgId": int(pdgIds[mother_idx]),
                            "mother_chain_indices": chain_indices,
                            "mother_chain_pdgIds": chain_pdgids,
                        })

            gen_objs.sort(key=lambda x: x["pt"], reverse=True)
            if len(gen_objs) < 2:
                append_gen_mother_records(gen_objs, set(), z, os.path.basename(iFile), iEv)
                continue

            # Find best OSSF gen pair near resonance
            best_pair = None
            best_delta = float("inf")
            for i, obj1 in enumerate(gen_objs):
                for obj2 in gen_objs[i+1:]:
                    if obj1["pdgId"] * obj2["pdgId"] >= 0:
                        continue
                    p41 = ROOT.TLorentzVector(); p41.SetPtEtaPhiM(obj1["pt"], obj1["eta"], obj1["phi"], obj1["mass"])
                    p42 = ROOT.TLorentzVector(); p42.SetPtEtaPhiM(obj2["pt"], obj2["eta"], obj2["phi"], obj2["mass"])
                    inv_mass = (p41 + p42).M()

                    if samp == "jpsi" and not (2.9 <= inv_mass <= 3.3):
                        continue
                    if samp == "dy" and not (60 <= inv_mass <= 120):
                        continue
                    if samp == "zprime" and not (2300 <= inv_mass <= 2700):
                        continue

                    delta = abs(inv_mass - resonance_mass)
                    if delta < best_delta:
                        best_delta = delta
                        best_pair = (obj1, obj2)

            if best_pair is None:
                append_gen_mother_records(gen_objs, set(), z, os.path.basename(iFile), iEv)
                continue

            best_indices = {best_pair[0]["index"], best_pair[1]["index"]}
            append_gen_mother_records(gen_objs, best_indices, z, os.path.basename(iFile), iEv)

            # Reco muons and matching
            if not hasattr(ev, "Muon_genPartIdx"):
                raise AttributeError("Muon_genPartIdx branch not found in the event tree.")
            muon_genPartIdx = getattr(ev, "Muon_genPartIdx")
            muon_genPartFlav = getattr(ev, "Muon_genPartFlav", None)
            muon_isTracker = getattr(ev, "Muon_isTracker", None)
            muon_isGlobal = getattr(ev, "Muon_isGlobal", None)
            muon_mediumId = getattr(ev, "Muon_mediumId", None)
            muon_pdgids = getattr(ev, "Muon_pdgId", None)
            muon_charge = getattr(ev, "Muon_charge", None)

            reco_muons = []
            for i_mu in range(ev.nMuon):
                gen_part_idx = -1
                if muon_genPartIdx is not None and i_mu < len(muon_genPartIdx):
                    gen_part_idx = int(muon_genPartIdx[i_mu])
                    if not (-1 <= gen_part_idx < nGenPart):
                        gen_part_idx = -1
                reco_muons.append({
                    "index": i_mu,
                    "pt": float(ev.Muon_pt[i_mu]),
                    "eta": float(ev.Muon_eta[i_mu]),
                    "phi": float(ev.Muon_phi[i_mu]),
                    "mass": float(ev.Muon_mass[i_mu]),
                    "is_tracker": bool(muon_isTracker[i_mu]) if muon_isTracker is not None else False,
                    "is_global": bool(muon_isGlobal[i_mu]) if muon_isGlobal is not None else False,
                    "mediumId": bool(muon_mediumId[i_mu]) if muon_mediumId is not None else False,
                    "pdgid": int(muon_pdgids[i_mu]) if muon_pdgids is not None else 0,
                    "charge": int(muon_charge[i_mu]) if muon_charge is not None else 0,
                    "genPartIdx": gen_part_idx,
                    "genPartFlav": int(muon_genPartFlav[i_mu]) if muon_genPartFlav is not None else None,
                })

            matched_reco_muons = []
            for gen_obj in (best_pair[0], best_pair[1]):
                gen_idx = int(gen_obj["index"])
                matching_muons = [mu for mu in reco_muons if mu["genPartIdx"] == gen_idx]
                matching_muons.sort(key=lambda m: m["pt"], reverse=True)
                matched_muon = matching_muons[0] if matching_muons else None

                tracker_pass = bool(matched_muon["is_tracker"]) if matched_muon else False
                global_pass  = bool(matched_muon["is_global"] and matched_muon["is_tracker"]) if matched_muon else False

                gen_status_records.append({
                    "pt": gen_obj["pt"],
                    "eta": gen_obj["eta"],
                    "phi": gen_obj["phi"],
                    "mass": gen_obj["mass"],
                    "from_target": True,
                    "assigned": matched_muon is not None,
                    "tracker_pass": tracker_pass,
                    "global_pass": global_pass,
                })

                if matched_muon:
                    matched_reco_muons.append(matched_muon)
                    muon_pdg = int(matched_muon.get("pdgid", 0))
                    if abs(muon_pdg) != 13:
                        reco_mismatch_records.append({
                            "file_index": z,
                            "file_name": os.path.basename(iFile),
                            "event_index": iEv,
                            "gen_index": gen_idx,
                            "gen_pdgId": int(gen_obj["pdgId"]),
                            "muon_index": matched_muon["index"],
                            "muon_pdgId": muon_pdg,
                            "muon_pt": float(matched_muon["pt"]),
                            "expected_abs_pdgId": 13,
                        })
                        lookup_key = (z, iEv, gen_idx)
                        rec_idx = record_lookup.get(lookup_key)
                        if rec_idx is not None:
                            gen_mother_records[rec_idx]["reco_mismatch_pdgId"] = muon_pdg

                    gen_pt = gen_obj["pt"]
                    if gen_pt > 0:
                        res = (matched_muon["pt"] - gen_pt) / gen_pt
                        reco_momentum_resolutions.append({
                            "gen_eta": float(gen_obj["eta"]),
                            "reco_eta": float(matched_muon["eta"]),
                            "gen_pt": float(gen_pt),
                            "reco_pt": float(matched_muon["pt"]),
                            "resolution": float(res),
                        })

            if len(matched_reco_muons) >= 2:
                matched_reco_muons.sort(key=lambda x: x["pt"], reverse=True)
                p41 = ROOT.TLorentzVector(); p41.SetPtEtaPhiM(
                    matched_reco_muons[0]["pt"], matched_reco_muons[0]["eta"], matched_reco_muons[0]["phi"], matched_reco_muons[0]["mass"]
                )
                p42 = ROOT.TLorentzVector(); p42.SetPtEtaPhiM(
                    matched_reco_muons[1]["pt"], matched_reco_muons[1]["eta"], matched_reco_muons[1]["phi"], matched_reco_muons[1]["mass"]
                )
                event_rows.append({
                    "mass": float((p41 + p42).M()),
                    "lead_pt": float(matched_reco_muons[0]["pt"]),
                    "lead_eta": float(matched_reco_muons[0]["eta"]),
                    "lead_phi": float(matched_reco_muons[0]["phi"]),
                    "sub_pt": float(matched_reco_muons[1]["pt"]),
                    "sub_eta": float(matched_reco_muons[1]["eta"]),
                    "sub_phi": float(matched_reco_muons[1]["phi"]),
                })

        tf.Close()

    print(f"Finished processing {nF} files.")
    print(f"Total selected gen {obj}s: {total_gen_selected}")
    print(f"Total selected gen {obj}s with {target_label} ancestor (|PDG|={target_pdg}): {total_gen_from_target}")

    # Build bin containers
    pt_bins = [(pt_edges[i], pt_edges[i+1]) for i in range(len(pt_edges)-1)]
    eta_bins = [(eta_edges[i], eta_edges[i+1]) for i in range(len(eta_edges)-1)]
    n_eta, n_pt = len(eta_bins), len(pt_bins)

    total_counts_gen = [[0 for _ in range(n_pt)] for _ in range(n_eta)]
    gen_to_tracker_counts = [[0 for _ in range(n_pt)] for _ in range(n_eta)]
    tracker_to_global_counts = [[0 for _ in range(n_pt)] for _ in range(n_eta)]
    tracker_to_global_total_counts = [[0 for _ in range(n_pt)] for _ in range(n_eta)]

    # Fill counts
    for rec in gen_status_records:
        if not rec.get("from_target", False):
            continue
        pt = rec["pt"]
        aeta = abs(rec["eta"])
        assigned = bool(rec.get("assigned", False))
        tracker_pass = bool(rec.get("tracker_pass", False))
        global_pass = bool(rec.get("global_pass", False))
        for ie, (eta_lo, eta_hi) in enumerate(eta_bins):
            if not (eta_lo <= aeta < eta_hi):
                continue
            for ip, (pt_lo, pt_hi) in enumerate(pt_bins):
                if not (pt_lo <= pt < pt_hi):
                    continue
                total_counts_gen[ie][ip] += 1
                if assigned and tracker_pass:
                    gen_to_tracker_counts[ie][ip] += 1
                    tracker_to_global_total_counts[ie][ip] += 1
                    if global_pass:
                        tracker_to_global_counts[ie][ip] += 1
                break

    # Compute efficiencies with Clopper-Pearson
    gen_to_tracker_eff = [[np.nan for _ in range(n_pt)] for _ in range(n_eta)]
    gen_to_tracker_err_lo = [[0.0 for _ in range(n_pt)] for _ in range(n_eta)]
    gen_to_tracker_err_up = [[0.0 for _ in range(n_pt)] for _ in range(n_eta)]

    tracker_to_global_eff = [[np.nan for _ in range(n_pt)] for _ in range(n_eta)]
    tracker_to_global_err_lo = [[0.0 for _ in range(n_pt)] for _ in range(n_eta)]
    tracker_to_global_err_up = [[0.0 for _ in range(n_pt)] for _ in range(n_eta)]

    for ie in range(n_eta):
        for ip in range(n_pt):
            k = gen_to_tracker_counts[ie][ip]
            n = total_counts_gen[ie][ip]
            if n > 0:
                eff = k / n
                low = ROOT.TEfficiency.ClopperPearson(n, k, alpha, False)
                up = ROOT.TEfficiency.ClopperPearson(n, k, alpha, True)
                gen_to_tracker_eff[ie][ip] = eff
                gen_to_tracker_err_lo[ie][ip] = eff - low
                gen_to_tracker_err_up[ie][ip] = up - eff

            k2 = tracker_to_global_counts[ie][ip]
            n2 = tracker_to_global_total_counts[ie][ip]
            if n2 > 0:
                eff2 = k2 / n2
                low2 = ROOT.TEfficiency.ClopperPearson(n2, k2, alpha, False)
                up2 = ROOT.TEfficiency.ClopperPearson(n2, k2, alpha, True)
                tracker_to_global_eff[ie][ip] = eff2
                tracker_to_global_err_lo[ie][ip] = eff2 - low2
                tracker_to_global_err_up[ie][ip] = up2 - eff2

    # CSV: efficiencies (write numbers directly; let NaN be "nan")
    with open(f"{samp}_{obj}_m2_efficiency.csv", "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow([
            "eta_lo","eta_hi","pt_lo","pt_hi",
            "gen_total","gen_matched_tracker","eff_gen_to_tracker","err_lo_gen_to_tracker","err_up_gen_to_tracker",
            "tracker_total","tracker_matched_global","eff_tracker_to_global","err_lo_tracker_to_global","err_up_tracker_to_global",
        ])
        for ie, (eta_lo, eta_hi) in enumerate(eta_bins):
            for ip, (pt_lo, pt_hi) in enumerate(pt_bins):
                writer.writerow([
                    eta_lo, eta_hi, pt_lo, pt_hi,
                    total_counts_gen[ie][ip],
                    gen_to_tracker_counts[ie][ip],
                    gen_to_tracker_eff[ie][ip],
                    gen_to_tracker_err_lo[ie][ip],
                    gen_to_tracker_err_up[ie][ip],
                    tracker_to_global_total_counts[ie][ip],
                    tracker_to_global_counts[ie][ip],
                    tracker_to_global_eff[ie][ip],
                    tracker_to_global_err_lo[ie][ip],
                    tracker_to_global_err_up[ie][ip],
                ])

    # CSV: gen ancestry
    with open(f"{samp}_{obj}_m2_gen_ancestry.csv", "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow([
            "file_index","file_name","event_index","gen_index","pdgId","pt","eta","phi","mass",
            "from_target","in_best_pair","mother_chain_indices","mother_chain_pdgIds",
            "mother_pdgId_if_not_target","reco_mismatch_pdgId",
        ])
        for rec in gen_mother_records:
            writer.writerow([
                rec["file_index"],
                rec["file_name"],
                rec["event_index"],
                rec["gen_index"],
                rec["pdgId"],
                rec["pt"],
                rec["eta"],
                rec["phi"],
                rec["mass"],
                int(rec["from_target"]),
                int(rec["in_best_pair"]),
                "->".join(map(str, rec["mother_chain_indices"])),
                "->".join(map(str, rec["mother_chain_pdgIds"])),
                "" if rec.get("mother_pdgId_non_target") is None else rec["mother_pdgId_non_target"],
                "" if rec.get("reco_mismatch_pdgId") is None else rec["reco_mismatch_pdgId"],
            ])

    # CSV: reco PDG mismatches (optional)
    if reco_mismatch_records:
        with open(f"{samp}_{obj}_m2_reco_mismatches.csv", "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow([
                "file_index","file_name","event_index","gen_index","gen_pdgId",
                "muon_index","muon_pdgId","muon_pt","expected_abs_pdgId",
            ])
            for rec in reco_mismatch_records:
                writer.writerow([
                    rec["file_index"], rec["file_name"], rec["event_index"], rec["gen_index"], rec["gen_pdgId"],
                    rec["muon_index"], rec["muon_pdgId"], rec["muon_pt"], rec["expected_abs_pdgId"],
                ])
        print(f"Saved reco mismatch CSV: {samp}_{obj}_m2_reco_mismatches.csv ({len(reco_mismatch_records)} entries)")
    else:
        print("No reco muon PDG ID mismatches found for matched leptons.")

    # CSV: per-event reco pair kinematics
    with open(f"{samp}_{obj}_m2_kins.csv", "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["mass","lead_pt","lead_eta","lead_phi","sub_pt","sub_eta","sub_phi"])
        for r in event_rows:
            writer.writerow([r['mass'], r['lead_pt'], r['lead_eta'], r['lead_phi'], r['sub_pt'], r['sub_eta'], r['sub_phi']])

    # CSV: per-match momentum resolution
    with open(f"{samp}_{obj}_m2_momentum_resolution_reco.csv", "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["gen_eta","reco_eta","gen_pt","reco_pt","resolution"])
        for entry in reco_momentum_resolutions:
            writer.writerow([entry['gen_eta'], entry['reco_eta'], entry['gen_pt'], entry['reco_pt'], entry['resolution']])

if __name__ == "__main__":
    main()

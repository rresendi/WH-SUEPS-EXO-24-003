import ROOT, os, sys, csv, math
ROOT.gROOT.SetBatch(True)

samp = sys.argv[1]
input_paths = sys.argv[2:]
require_medium = False

# Collect root files
files = []
for p in input_paths:
    if os.path.isdir(p):
        for fn in os.listdir(p):
            if fn.endswith(".root"):
                files.append(os.path.join(p, fn))
    elif p.endswith(".root"):
        files.append(p)
if not files:
    print("No ROOT files found.")
    sys.exit(1)

# Binning and sample setup
eta_edges = [0.0, 0.9, 1.2, 2.1, 2.4]
if samp == "jpsi_b":
    eta_edges = [0.0, 1.2]
    pt_edges = [1.0, 5.0, 10.0]
    res_mass = 3.0969
    mwin = (2.9, 3.3)
elif samp == "jpsi_e":
    eta_edges = [1.2, 2.4]
    pt_edges = [1.0, 3.0, 10.0]
    res_mass = 3.0969
    mwin = (2.9, 3.3)
elif samp == "dy":
    pt_edges = [10.0, 15.0, 100.0, 1000.0]
    res_mass = 91.1876
    mwin = (60.0, 120.0)
else:  # zprime
    pt_edges = [1000.0, 3000.0]
    res_mass = 2500.0
    mwin = (2300.0, 2700.0)

pt_bins = [(pt_edges[i], pt_edges[i+1]) for i in range(len(pt_edges)-1)]
eta_bins = [(eta_edges[i], eta_edges[i+1]) for i in range(len(eta_edges)-1)]
n_eta, n_pt = len(eta_bins), len(pt_bins)

# Counters
gen_totals   = [[0 for _ in range(n_pt)] for __ in range(n_eta)]
gen_matched  = [[0 for _ in range(n_pt)] for __ in range(n_eta)]

# Per-object and per-event outputs
res_rows = []
event_rows = []
denominator_rows = []

# Clopper–Pearson uncertainty setup
CL = 0.683
ALPHA = 1.0 - CL
def cp_bounds(n, k):
    if n <= 0:
        return (float("nan"), 0.0, 0.0)
    eff = k / n
    low = ROOT.TEfficiency.ClopperPearson(n, k, ALPHA, False)
    up  = ROOT.TEfficiency.ClopperPearson(n, k, ALPHA, True)
    return (eff, eff - low, up - eff)

# Best OSSF gen pair
def pick_best_pair(gen_list):
    best = None
    best_delta = 1e99
    for i in range(len(gen_list)):
        a = gen_list[i]
        for j in range(i+1, len(gen_list)):
            b = gen_list[j]
            if a["pdgId"] * b["pdgId"] >= 0:
                continue
            p41 = ROOT.TLorentzVector(); p41.SetPtEtaPhiM(a["pt"], a["eta"], a["phi"], a["mass"])
            p42 = ROOT.TLorentzVector(); p42.SetPtEtaPhiM(b["pt"], b["eta"], b["phi"], b["mass"])
            m = (p41 + p42).M()
            if m < mwin[0] or m > mwin[1]:
                continue
            d = abs(m - res_mass)
            if d < best_delta:
                best_delta = d
                best = (a, b)
    return best

# Loop files
for ix, path in enumerate(files):
    # print(f"[{ix+1}/{len(files)}] {path}")
    tf = ROOT.TFile.Open(path, "READ")
    if not tf or tf.IsZombie():
        continue
    evs = tf.Get("Events")
    if not evs:
        tf.Close()
        continue

    for iev, ev in enumerate(evs):
        nG = getattr(ev, "nGenPart", 0)
        pdg = ev.GenPart_pdgId
        pt  = ev.GenPart_pt
        eta = ev.GenPart_eta
        phi = ev.GenPart_phi
        sta = ev.GenPart_status
        mas = ev.GenPart_mass

        # Select gen objs
        gen_list = []
        pt_min, pt_max = pt_edges[0], pt_edges[-1]
        for i in range(nG):
            if abs(pdg[i]) != 13: continue
            if abs(eta[i]) >= 2.4: continue
            if int(sta[i]) != 1: continue
            if not (pt_min <= pt[i] < pt_max): continue
            gen_list.append({
                "index": i, "pt": float(pt[i]), "eta": float(eta[i]),
                "phi": float(phi[i]), "mass": float(mas[i]), "pdgId": int(pdg[i]),
            })

        if len(gen_list) < 2: continue
        gen_list.sort(key=lambda x: x["pt"], reverse=True)
        pair = pick_best_pair(gen_list)
        if not pair: continue

        # Select reco objs
        nMu = getattr(ev, "nMuon", 0)
        gidx = getattr(ev, "Muon_genPartIdx", None)
        mid  = getattr(ev, "Muon_mediumId", None)

        reco = []
        for i in range(nMu):
            gi = int(gidx[i]) if (gidx is not None and i < len(gidx)) else -1
            reco.append({
                "index": i,
                "pt": float(ev.Muon_pt[i]),
                "eta": float(ev.Muon_eta[i]),
                "phi": float(ev.Muon_phi[i]),
                "mass": float(ev.Muon_mass[i]),
                "genPartIdx": gi,
                "mediumId": bool(mid[i]) if (mid is not None and i < len(mid)) else False,
            })

        # Count & match
        matched_reco = []
        seen_gen = set() 
        used_reco = set()
        for g in (pair[0], pair[1]):
            gi = g["index"]
            if gi in seen_gen:
                continue
            seen_gen.add(gi)
            cand = [m for m in reco if m["genPartIdx"] == gi and (not require_medium or m["mediumId"]) and (m["index"] not in used_reco)]
            cand.sort(key=lambda x: x["pt"], reverse=True)
            mbest = cand[0] if cand else None

            geta = g["eta"]; aeta = abs(geta); gpt = g["pt"]
            for ie, (elo, ehi) in enumerate(eta_bins):
                if not (elo <= aeta < ehi): continue
                for ip, (plo, phi_) in enumerate(pt_bins):
                    if not (plo <= gpt < phi_): continue
                    gen_totals[ie][ip] += 1
                    denominator_rows.append((gpt, geta))
                    if mbest:
                        gen_matched[ie][ip] += 1
                        matched_reco.append(mbest)
                        used_reco.add(mbest["index"])
                        # compute deltaR and relative pT
                        dphi = abs(mbest["phi"] - g["phi"])
                        if dphi > math.pi: dphi = 2*math.pi - dphi
                        dR = math.hypot(mbest["eta"] - g["eta"], dphi)
                        rel_pt = (mbest["pt"] - gpt)/gpt if gpt > 0 else float("nan")
                        if gpt > 0:
                            res = (mbest["pt"] - gpt) / gpt
                            res_rows.append([g["eta"], mbest["eta"], gpt, mbest["pt"], res])
                    break
                break

        # Reco mass
        if len(matched_reco) >= 2:
            matched_reco.sort(key=lambda x: x["pt"], reverse=True)
            p41 = ROOT.TLorentzVector(); p41.SetPtEtaPhiM(matched_reco[0]["pt"], matched_reco[0]["eta"], matched_reco[0]["phi"], matched_reco[0]["mass"])
            p42 = ROOT.TLorentzVector(); p42.SetPtEtaPhiM(matched_reco[1]["pt"], matched_reco[1]["eta"], matched_reco[1]["phi"], matched_reco[1]["mass"])
            m = (p41 + p42).M()
            event_rows.append([m,
                               matched_reco[0]["pt"], matched_reco[0]["eta"], matched_reco[0]["phi"],
                               matched_reco[1]["pt"], matched_reco[1]["eta"], matched_reco[1]["phi"]])

    tf.Close()


# ---- Write CSVs ----
if require_medium == True:
    prefix = f"{samp}_muon_id"
else:
    prefix = f"{samp}_muon_noid"

with open(f"{prefix}_efficiency.csv", "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["eta_lo","eta_hi","pt_lo","pt_hi","gen_total","gen_matched","eff","err_lo","err_up"])
    for ie, (elo, ehi) in enumerate(eta_bins):
        for ip, (plo, phi_) in enumerate(pt_bins):
            n = gen_totals[ie][ip]; k = gen_matched[ie][ip]
            eff, el, eu = cp_bounds(n, k)
            w.writerow([elo, ehi, plo, phi_, n, k,
                        f"{eff:.6f}" if n>0 else "",
                        f"{el:.6f}" if n>0 else "",
                        f"{eu:.6f}" if n>0 else ""])

with open(f"{prefix}_momentum_resolution.csv", "w", newline="") as f:
    w = csv.writer(f); w.writerow(["gen_eta","reco_eta","gen_pt","reco_pt","resolution"])
    for r in res_rows:
        w.writerow(r)

with open(f"{prefix}_kins.csv", "w", newline="") as f:
    w = csv.writer(f); w.writerow(["mass","lead_pt","lead_eta","lead_phi","sub_pt","sub_eta","sub_phi"])
    for r in event_rows:
        w.writerow(r)

with open(f"{prefix}_gen_denominator.csv", "w", newline="") as f:
    w = csv.writer(f); w.writerow(["gen_pt","gen_eta"])
    for rec in denominator_rows:
        w.writerow(rec)

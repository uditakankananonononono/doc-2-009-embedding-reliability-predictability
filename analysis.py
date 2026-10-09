"""DOC-2-009 R3 frozen analysis (PROTOCOL.md lock-1). Run once. Needs per_assay.csv, exposures.tsv."""
import json, numpy as np, pandas as pd
rng = np.random.default_rng(12345)
a = pd.read_csv("per_assay.csv"); e = pd.read_csv("exposures.tsv", sep="\t").drop(columns=["seq_len"])
d = a.merge(e, on="DMS_id"); d = d[d.rho.notna() & (d.n_scored >= 50)].reset_index(drop=True)
d["tsub"] = d.DMS_id.str.contains("Tsuboyama_2023").astype(float); d["loglen"] = np.log(d.seq_len)
d["ann_missing"] = d.annotation_score.isna().astype(float); d["ann"] = d.annotation_score.fillna(0.0); d["lp"] = np.log1p(d.n_pubmed.fillna(0))
ov = set(d[d.year <= 2021].UniProt_ID) & set(d[d.year >= 2022].UniProt_ID)
disc = d[d.year <= 2021].reset_index(drop=True); ext = d[(d.year >= 2022) & (~d.UniProt_ID.isin(ov))].reset_index(drop=True)
def X(df, kind):
    c = []
    if kind in ("B1", "M"): c += [(df.taxon == t).astype(float).values for t in ["Eukaryote", "Prokaryote", "Virus"]] + [df.tsub.values]
    if kind == "M": c += [df.loglen.values, df.wtll.values, df.ann.values, df.ann_missing.values, df.lp.values]
    return np.column_stack(c) if c else np.zeros((len(df), 0))
def fit_pred(Xtr, ytr, Xte, lam=1.0):
    if Xtr.shape[1] == 0: return np.full(len(Xte), ytr.mean())
    mu, sd = Xtr.mean(0), Xtr.std(0) + 1e-9; A = (Xtr - mu) / sd; B = (Xte - mu) / sd
    w = np.linalg.solve(A.T @ A + lam * np.eye(A.shape[1]), A.T @ (ytr - ytr.mean())); return B @ w + ytr.mean()
def lopo(df, kind):
    pred = np.zeros(len(df)); ref = np.zeros(len(df)); Xa = X(df, kind)
    for p in df.UniProt_ID.unique():
        te = (df.UniProt_ID == p).values; tr = ~te
        pred[te] = fit_pred(Xa[tr], df.rho.values[tr], Xa[te]); ref[te] = df.rho.values[tr].mean()
    y = df.rho.values; return 1 - ((y - pred) ** 2).sum() / ((y - ref) ** 2).sum()
def ext_r2(tr, te, kind, shuffle_map=None):
    Xt, Xe = X(tr, kind), X(te, kind)
    if shuffle_map is not None: Xt = shuffle_map(tr, Xt)
    p = fit_pred(Xt, tr.rho.values, Xe); y = te.rho.values
    return 1 - ((y - p) ** 2).sum() / ((y - tr.rho.mean()) ** 2).sum()
res = {"n_disc": len(disc), "n_ext": len(ext), "overlap_removed": sorted(ov)}
r = {k: lopo(disc, k) for k in ["B0", "B1", "M"]}; res["G1"] = dict(lopo_r2=r, delta=r["M"] - r["B1"], pass_=bool(r["M"] - r["B1"] >= 0.05))
eb1, em = ext_r2(disc, ext, "B1"), ext_r2(disc, ext, "M"); delta = em - eb1
prots = ext.UniProt_ID.unique(); grp = {p: ext[ext.UniProt_ID == p] for p in prots}; bs = []
for _ in range(10000):
    s = pd.concat([grp[p] for p in rng.choice(prots, len(prots))]).reset_index(drop=True)
    y = s.rho.values; Xs1, Xsm = X(s, "B1"), X(s, "M"); Xd1, Xdm = X(disc, "B1"), X(disc, "M")
    p1 = fit_pred(Xd1, disc.rho.values, Xs1); pm = fit_pred(Xdm, disc.rho.values, Xsm); den = ((y - disc.rho.mean()) ** 2).sum()
    if den > 0: bs.append((((y - p1) ** 2).sum() - ((y - pm) ** 2).sum()) / den)
ci = list(np.percentile(bs, [2.5, 97.5]))
res["G2"] = dict(ext_r2_M=em, ext_r2_B1=eb1, delta=delta, ci=ci, pass_=bool(em > 0 and delta >= 0.05 and ci[0] > 0))
base_cols = ["loglen", "wtll", "ann", "ann_missing", "lp"]; pu = disc.UniProt_ID.unique(); perm = []
for _ in range(2000):
    sh = pd.Series(rng.permutation(len(pu)), index=pu); keyed = {u: i for i, u in enumerate(pu)}
    dd = disc.copy(); idx = dd.UniProt_ID.map(lambda u: pu[sh[u]]); first = disc.drop_duplicates("UniProt_ID").set_index("UniProt_ID")[base_cols]
    for c in base_cols: dd[c] = idx.map(first[c]).values
    perm.append(ext_r2(dd, ext, "M") - eb1)
pv = (1 + sum(x >= delta for x in perm)) / (1 + len(perm)); res["G3"] = dict(perm_p=pv, pass_=bool(pv < 0.05))
res["LABEL"] = "POSITIVE" if all(res[g]["pass_"] for g in ["G1", "G2", "G3"]) else "HONEST NEGATIVE"
print("RESULT_JSON", json.dumps(res, default=float)); open("results.json", "w").write(json.dumps(res, default=float, indent=1))

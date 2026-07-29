import numpy as np


def ref_distances(store, cohort_pids, refs, arm, case=1, control=0):
    coh = store.meta.loc[cohort_pids]
    centroid = {}
    for lab in (case, control):
        pids = coh.index[coh.label == lab].tolist()
        if pids:
            centroid[lab] = store.sigs(pids).mean(0)
    sig = {r: store.sigs([r])[0] for r in refs}
    lbl = {r: int(store.meta.loc[r, "label"]) for r in refs}
    rows = []
    for r in refs:
        target = lbl[r] if arm == "cow" else (control if lbl[r] == case else case)
        d_cohort = float(np.linalg.norm(sig[r] - centroid[target])) if target in centroid else float("nan")
        same = [o for o in refs if o != r and lbl[o] == lbl[r]]
        d_refs = float(min(np.linalg.norm(sig[r] - sig[o]) for o in same)) if same else float("nan")
        rows.append({"patient": r, "study": str(store.meta.loc[r, "study"]),
                     "label": lbl[r], "d_cohort": d_cohort, "d_refs": d_refs})
    return rows

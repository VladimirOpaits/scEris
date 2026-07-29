import argparse
import glob
import os
import numpy as np
import scanpy as sc


def expand(src):
    path, _, gcol = src.partition(":")
    files = sorted(glob.glob(os.path.join(path, "*.h5ad"))) if os.path.isdir(path) else [path]
    return [(f, gcol or None) for f in files]


def stage(path, gene_col, out, celltype_col, donor_col, per_donor, max_per_type, seed):
    a = sc.read_h5ad(path, backed="r")
    genes = (a.var[gene_col] if gene_col else a.var.index).astype(str).values
    rng = np.random.default_rng(seed)
    cts = a.obs[celltype_col].to_numpy()
    donors = a.obs[donor_col].to_numpy()
    keep = []
    for ct in np.unique(cts):
        pool = np.where(cts == ct)[0]
        d = donors[pool]
        chosen = np.concatenate([
            rng.choice(pool[d == dv], min((d == dv).sum(), per_donor), replace=False)
            for dv in np.unique(d)
        ])
        if len(chosen) > max_per_type:
            chosen = rng.choice(chosen, max_per_type, replace=False)
        keep.append(chosen)
    keep = np.sort(np.concatenate(keep))
    sub = a[keep].to_memory()
    sub.var_names = genes
    name = os.path.splitext(os.path.basename(path))[0]
    sub.write_h5ad(os.path.join(out, f"{name}.h5ad"))
    print(f"{name}: {sub.n_obs} cells staged", flush=True)


def main():
    ap = argparse.ArgumentParser(prog="sceris-stage")
    ap.add_argument("--src", action="append", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--celltype-col", default="cell_type")
    ap.add_argument("--donor-col", default="donor_id")
    ap.add_argument("--per-donor", type=int, default=50)
    ap.add_argument("--max-per-type", type=int, default=5000)
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)
    for src in args.src:
        for path, gcol in expand(src):
            stage(path, gcol, args.out, args.celltype_col, args.donor_col, args.per_donor, args.max_per_type, args.seed)


if __name__ == "__main__":
    main()

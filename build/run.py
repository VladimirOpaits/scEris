import argparse
import glob
import os
import anndata as ad
import scanpy as sc
from basis import Basis


def load_dir(path, gene_col=None):
    files = sorted(glob.glob(os.path.join(path, "*.h5ad")))
    if not files:
        raise FileNotFoundError(f"no .h5ad in {path}")
    parts = []
    for f in files:
        a = sc.read_h5ad(f)
        if gene_col:
            a.var_names = a.var[gene_col].astype(str).values
        parts.append(a)
    a = parts[0] if len(parts) == 1 else ad.concat(parts, join="inner")
    a.obs_names_make_unique()
    return a


def select_hvg(adata, n):
    t = adata.copy()
    sc.pp.normalize_total(t, target_sum=1e4)
    sc.pp.log1p(t)
    sc.pp.highly_variable_genes(t, n_top_genes=n)
    return adata[:, t.var.highly_variable].copy()


def _safe(s):
    return "".join(c if c.isalnum() else "_" for c in str(s))


def build_one(adata, args):
    a = select_hvg(adata, args.n_hvg) if args.n_hvg else adata
    return Basis(n_pca=args.n_pca, chunk=args.chunk).fit(a)


def main():
    ap = argparse.ArgumentParser(prog="sceris-build")
    ap.add_argument("--dir", required=True)
    ap.add_argument("--out", default="basis.npz")
    ap.add_argument("--gene-col", default=None)
    ap.add_argument("--by-celltype", default=None)
    ap.add_argument("--compartments", action="store_true")
    ap.add_argument("--cl-obo", default=None)
    ap.add_argument("--celltype-col", default="cell_type")
    ap.add_argument("--cl-col", default="cell_type_ontology_term_id")
    ap.add_argument("--min-cells", type=int, default=500)
    ap.add_argument("--n-hvg", type=int, default=None)
    ap.add_argument("--n-pca", type=int, default=50)
    ap.add_argument("--chunk", type=int, default=20000)
    args = ap.parse_args()

    adata = load_dir(args.dir, args.gene_col)
    print(f"loaded {adata.n_obs} cells x {adata.n_vars} genes", flush=True)

    group_col = args.by_celltype
    if args.compartments:
        from compartments import load_graph, build_mapper, CL_URL
        mapper = build_mapper(load_graph(args.cl_obo or CL_URL))
        cids = adata.obs[args.cl_col].astype(str).values if args.cl_col in adata.obs else [None] * adata.n_obs
        cts = adata.obs[args.celltype_col].astype(str).values
        adata.obs["compartment"] = [mapper(ct, ci) for ct, ci in zip(cts, cids)]
        group_col = "compartment"
        print("compartments:", dict(adata.obs["compartment"].value_counts()), flush=True)

    if group_col:
        os.makedirs(args.out, exist_ok=True)
        for ct in sorted(adata.obs[group_col].astype(str).unique()):
            sub = adata[adata.obs[group_col].astype(str) == ct].copy()
            if sub.n_obs < args.min_cells:
                print(f"skip {ct}: {sub.n_obs} cells", flush=True)
                continue
            build_one(sub, args).save(os.path.join(args.out, _safe(ct) + ".npz"))
            print(f"{ct}: {sub.n_obs} cells -> saved", flush=True)
    else:
        build_one(adata, args).save(args.out)
        print(f"basis saved -> {args.out}", flush=True)


if __name__ == "__main__":
    main()

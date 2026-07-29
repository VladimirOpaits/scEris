import numpy as np
import scipy.sparse as sp


def load_basis(path):
    d = np.load(path)
    return {k: d[k] for k in d.files}


def _normalize(counts, target):
    X = counts.tocsr().astype(np.float32)
    rs = np.asarray(X.sum(1)).ravel()
    rs[rs == 0] = 1.0
    X.data *= np.repeat((target / rs).astype(np.float32), np.diff(X.indptr))
    X.data = np.log1p(X.data)
    return X


def apply_basis(basis, counts, var_names, chunk=20000):
    genes = [str(g) for g in basis["genes"]]
    pos = {g: i for i, g in enumerate(map(str, var_names))}
    idx = [pos[g] for g in genes]
    X = _normalize(counts[:, idx], float(basis["target"]))
    mean, scale, clip = basis["mean"], basis["scale"], float(basis["clip"])
    pca_mean, comp = basis["pca_mean"], basis["components"]
    out = []
    for i in range(0, X.shape[0], chunk):
        z = X[i:i + chunk].toarray()
        z = (z - mean) / scale
        np.clip(z, -clip, clip, out=z)
        out.append((z - pca_mean) @ comp.T)
    return np.vstack(out)

import numpy as np


def load_basis(path):
    d = np.load(path)
    return {k: d[k] for k in d.files}


def apply_basis(basis, counts, var_names, chunk=20000):
    genes = [str(g) for g in basis["genes"]]
    pos = {g: i for i, g in enumerate(map(str, var_names))}
    cols = [pos[g] for g in genes if g in pos]
    slots = [j for j, g in enumerate(genes) if g in pos]
    mean, scale, clip = basis["mean"], basis["scale"], float(basis["clip"])
    pca_mean, comp, target = basis["pca_mean"], basis["components"], float(basis["target"])
    counts = counts.tocsr()
    out = []
    for i in range(0, counts.shape[0], chunk):
        block = counts[i:i + chunk][:, cols].toarray().astype(np.float32)
        z = np.zeros((block.shape[0], len(genes)), np.float32)
        z[:, slots] = block
        rs = z.sum(1)
        rs[rs == 0] = 1.0
        z = np.log1p(z / rs[:, None] * target)
        z = (z - mean) / scale
        np.clip(z, -clip, clip, out=z)
        out.append((z - pca_mean) @ comp.T)
    return np.vstack(out)

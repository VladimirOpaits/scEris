import scanpy as sc
import scipy.sparse as sp
import numpy as np
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import IncrementalPCA
from tqdm import tqdm

def _dense(x):
    return x.toarray() if sp.issparse(x) else np.asarray(x)

class Basis:
    def __init__(self, n_pca = 50, clip = 10.0, target = 1e4, chunk = 20000):
        self.n_pca = n_pca
        self.clip = clip
        self.target = target
        self.chunk = chunk
        self.scaler = None
        self.ipca = None
        self.genes = None
        
    def _norm(self, adata):
        a = adata.copy()
        sc.pp.normalize_total(a, target_sum = self.target)       # считаем пропорции генов и умножаем на 1e4
        sc.pp.log1p(a)                                           # добавляет 1 и сжимает ближе к нормальному распределению
        return a.X

    def _scale(self, z):
        z = (z - self.scaler.mean_) / self.scaler.scale_
        np.clip(z, -self.clip, self.clip, out = z)
        return z
    
    def fit(self, adata, rows = None):
        self.genes = list(adata.var_names)
        X = self._norm(adata)
        if rows is not None:
            X = X[rows]
        n = X.shape[0]                                              # число клеток
        self.scaler = StandardScaler()
        for i in tqdm(range(0, n, self.chunk), desc = "scaler"):                           # потоковый расчет стандартизации
            self.scaler.partial_fit(_dense(X[i:i + self.chunk]))    # узнаем mean/std
        self.ipca = IncrementalPCA(n_components=self.n_pca)
        for i in tqdm(range(0, n, self.chunk), desc = "pca"):
            z = self._scale(_dense(X[i:i + self.chunk]))
            if z.shape[0] >= self.n_pca:                            # проверка, потому что нельзя получить 50 компонент из меньше чем 50 точек
                self.ipca.partial_fit(z)
        return self

    def transform(self, adata):
        X = self._norm(adata[:, self.genes])
        out = []
        for i in range(0, X.shape[0], self.chunk):
            out.append(self.ipca.transform(self._scale(_dense(X[i:i + self.chunk]))))
        return np.vstack(out)

    def save(self, path):
        with open(path, "wb") as f:
            np.savez(f, genes=np.array(self.genes), mean=self.scaler.mean_,
                     scale=self.scaler.scale_, clip=self.clip, target=self.target,
                     pca_mean=self.ipca.mean_, components=self.ipca.components_)
    
    
    
    
    
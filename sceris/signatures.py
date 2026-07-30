import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist
from sklearn.kernel_approximation import RBFSampler

def _median_sigma(X, n = 5000, seed = 0):
    rng = np.random.default_rng(seed)
    idx = rng.choice(len(X), min(n, len(X)), replace = False)  # выбираем случайные 5000 точек, что бы не считать медиану по всему что есть
    d = pdist(X[idx])    # pdist вычисляет попарные расстояния между всеми возможными точками, возвращает массив
    return float(np.median(d[d > 0]))     # проверка на ненужные нам нули

class Signatures:
    def __init__(self, rff_dim = 1024, seed = 0):
        self.rff_dim = rff_dim
        self.seed = seed
        
    def fit(self, coords):
        sigma = _median_sigma(coords, seed = self.seed)
        self.rbf = RBFSampler(gamma = 1.0 / (2 * sigma ** 2), n_components = self.rff_dim, random_state = self.seed)
        self.rbf.fit(coords)    # генерируем рандомную матрицу один раз
        return self             # цепочный вызов
    
    def patient(self, coords, patient_ids):
        phi = self.rbf.transform(coords)
        df = pd.DataFrame(phi)
        df["patient_id"] = np.asarray(patient_ids)
        return df.groupby("patient_id").mean()
    
    
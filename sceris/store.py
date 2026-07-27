import numpy as np
import pandas as pd

class SignatureStore:
    def __init__(self, sig: np.ndarray, meta: pd.DataFrame):
        assert len(sig) == len(meta)
        self.sig = np.ascontiguousarray(sig, dtype = np.float32)
        self.meta = meta
        self._row = {pid: i for i, pid in enumerate(meta.index)}
        
    def sigs(self, pids):
        return self.sig[[self._row[p] for p in pids]]
    
    def __len__(self):
        return len(self.meta)
    
    @classmethod         # фабрика 
    def from_parquet(cls, path, label_cols = "sample_type", pos_label = "tumor", study_col = "study", assay_col = "assay", id_col = "sample_id", sig_prefix = "s"):
        df = pd.read_parquet(path)
        sig_cols = [c for c in df.columns if c.startswith(sig_prefix) and c[len(sig_prefix):].isdigit()]  # все что идет после префикса должно быть числом
        sig = df[sig_cols].to_numpy(dtype = np.float32)        # лежит нампай матрица
        raw = df[label_cols].to_numpy()                        # сырая метка (tumor/normal/подтипы) — не теряем на входе
        label = (raw == pos_label).astype(int) if pos_label is not None else raw  # бинарный вью только если задан pos_label; иначе метка как есть (мультикласс)
        meta = pd.DataFrame({
            "label": label,
            "label_raw": raw,
            "study": df[study_col].to_numpy(),
            "assay": df[assay_col].to_numpy(),
            "n_cells": df["n_cells"].to_numpy(),
        }, index=pd.Index(df[id_col], name="patient_id"))
        return cls(sig, meta)
        
    
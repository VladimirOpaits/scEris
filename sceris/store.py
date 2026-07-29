import numpy as np
import pandas as pd
from pathlib import Path

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
        if "label" in df.columns and "patient_id" in df.columns:
            raw = label = df["label"].to_numpy()
            ids = df["patient_id"]
        else:
            raw = df[label_cols].to_numpy()
            label = (raw == pos_label).astype(int) if pos_label is not None else raw  # бинарный вид только если задан pos_label, иначе метка как есть
            ids = df[id_col]
        meta = pd.DataFrame({
            "label": label,
            "label_raw": raw,
            "study": df[study_col].to_numpy(),
            "assay": df[assay_col].to_numpy() if assay_col in df.columns else "?",
            "n_cells": df["n_cells"].to_numpy() if "n_cells" in df.columns else 0,
        }, index=pd.Index(ids, name="patient_id"))
        return cls(sig, meta)

    @classmethod
    def from_path(cls, path, **kw):
        p = Path(path)
        files = sorted(p.glob("*.parquet")) if p.is_dir() else [p]
        if not files:
            raise FileNotFoundError(f"no parquet under {path}")
        store = cls.from_parquet(str(files[0]), **kw)
        for f in files[1:]:
            store = store.extended(cls.from_parquet(str(f), **kw))
        return store

    def extended(self, other):             # собирает в один обьект, что бы можно было передавать все что надо вместе
        return SignatureStore(np.vstack([self.sig, other.sig]),
                              pd.concat([self.meta, other.meta]))

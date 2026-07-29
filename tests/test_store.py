import numpy as np
import pandas as pd
import pytest
from sceris.store import SignatureStore


def _mk(path, ids, d=4, seed=0):
    rng = np.random.default_rng(seed)
    df = pd.DataFrame(rng.standard_normal((len(ids), d)), columns=[f"s{i}" for i in range(d)])
    df.insert(0, "patient_id", ids)
    df.insert(1, "label", [i % 2 for i in range(len(ids))])
    df.insert(2, "study", [f"S{i % 2}" for i in range(len(ids))])
    df.to_parquet(path)
    return path


def test_canonical_schema(tmp_path):
    s = SignatureStore.from_parquet(_mk(tmp_path / "a.parquet", ["p0", "p1", "p2", "p3"]))
    assert list(s.meta.columns) == ["label", "label_raw", "study", "assay", "n_cells"]
    assert set(s.meta.label.unique()) == {0, 1}
    assert (s.meta.assay == "?").all() and (s.meta.n_cells == 0).all()
    assert s.sig.shape == (4, 4)


def test_from_path_file(tmp_path):
    p = _mk(tmp_path / "a.parquet", ["p0", "p1"])
    assert len(SignatureStore.from_path(p)) == 2


def test_from_path_dir(tmp_path):
    _mk(tmp_path / "a.parquet", ["p0", "p1", "p2"])
    _mk(tmp_path / "b.parquet", ["q0", "q1"])
    assert len(SignatureStore.from_path(tmp_path)) == 5


def test_from_path_empty(tmp_path):
    with pytest.raises(FileNotFoundError):
        SignatureStore.from_path(tmp_path)


def test_extended(tmp_path):
    a = SignatureStore.from_parquet(_mk(tmp_path / "a.parquet", ["p0", "p1"], seed=1))
    b = SignatureStore.from_parquet(_mk(tmp_path / "b.parquet", ["q0", "q1", "q2"], seed=2))
    c = a.extended(b)
    assert len(c) == 5
    assert c.sig.shape == (5, 4)
    assert np.allclose(c.sigs(["q0"]), b.sigs(["q0"]))
    assert a.sig.shape == (2, 4)

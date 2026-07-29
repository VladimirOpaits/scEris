import numpy as np
import pandas as pd
import pytest
from sceris.store import SignatureStore
from sceris.stats import ref_distances


@pytest.fixture
def store():
    sig = np.array([[0., 0], [0, 0], [10, 0], [10, 0], [1, 0], [2, 0], [11, 0]], dtype=np.float32)
    meta = pd.DataFrame({
        "label": [1, 1, 0, 0, 1, 1, 0],
        "study": ["A", "A", "B", "B", "R", "R", "R"],
    }, index=pd.Index(["cA", "cB", "cC", "cD", "r1", "r2", "r3"], name="patient_id"))
    return SignatureStore(sig, meta)


COH = ["cA", "cB", "cC", "cD"]


def test_pig_d_cohort_is_opposite_class(store):
    rows = {r["patient"]: r for r in ref_distances(store, COH, ["r1", "r2"], "pig")}
    assert rows["r1"]["d_cohort"] == pytest.approx(9.0)
    assert rows["r2"]["d_cohort"] == pytest.approx(8.0)


def test_cow_d_cohort_is_same_class(store):
    rows = {r["patient"]: r for r in ref_distances(store, COH, ["r1", "r2"], "cow")}
    assert rows["r1"]["d_cohort"] == pytest.approx(1.0)
    assert rows["r2"]["d_cohort"] == pytest.approx(2.0)


def test_d_refs_same_class_nearest(store):
    rows = {r["patient"]: r for r in ref_distances(store, COH, ["r1", "r2"], "cow")}
    assert rows["r1"]["d_refs"] == pytest.approx(1.0)
    assert rows["r2"]["d_refs"] == pytest.approx(1.0)


def test_d_refs_nan_when_alone_in_class(store):
    rows = {r["patient"]: r for r in ref_distances(store, COH, ["r1", "r3"], "cow")}
    assert np.isnan(rows["r1"]["d_refs"])
    assert np.isnan(rows["r3"]["d_refs"])

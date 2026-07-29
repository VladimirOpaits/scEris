import os
import numpy as np
import pandas as pd
import pytest
from sceris import project
from sceris.source import Source
from sceris.dataset import load_cohort, spec_for


def _mk(path, ids=("p0", "p1", "p2"), d=4):
    rng = np.random.default_rng(0)
    df = pd.DataFrame(rng.standard_normal((len(ids), d)), columns=[f"s{i}" for i in range(d)])
    df.insert(0, "patient_id", list(ids))
    df.insert(1, "label", [i % 2 for i in range(len(ids))])
    df.insert(2, "study", [f"S{i % 2}" for i in range(len(ids))])
    df.to_parquet(path)
    return path


def test_source_local_absolute():
    spec = Source.local("data/x.parquet").spec
    assert spec["kind"] == "local"
    assert os.path.isabs(spec["path"])


def test_source_load(tmp_path):
    store = Source.local(_mk(tmp_path / "c.parquet")).load()
    assert len(store) == 3


def test_source_unknown_kind():
    with pytest.raises(ValueError):
        Source({"kind": "nope"}).load()


def test_source_census_not_implemented():
    with pytest.raises(NotImplementedError):
        Source({"kind": "census"}).load()


def test_dataset_parquet(tmp_path):
    store = load_cohort(spec_for(_mk(tmp_path / "d.parquet")))
    assert len(store) == 3


def test_dataset_cells_needs_source():
    with pytest.raises(ValueError):
        load_cohort({"kind": "cells", "path": "x"})


def test_dataset_unknown():
    with pytest.raises(ValueError):
        load_cohort({"kind": "nope"})


def test_spec_for_absolute():
    spec = spec_for("data/x.parquet")
    assert spec["kind"] == "parquet" and os.path.isabs(spec["path"])


def test_project_roundtrip(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    assert project.load() == {}
    project.set("source", {"kind": "local", "path": "/abs/x"})
    project.set("dataset", {"kind": "parquet", "path": "/abs/y"})
    assert project.get("source") == {"kind": "local", "path": "/abs/x"}
    assert set(project.load()) == {"source", "dataset"}

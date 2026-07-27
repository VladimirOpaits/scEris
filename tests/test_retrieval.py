from sceris import SignatureStore, diverse_refs, matched_refs, config
from sceris.retrieval import diverse_refs_paired, representatives
import pytest


@pytest.fixture(scope = "session")
def store():
    return SignatureStore.from_parquet(path=config.CRC_SIG)

@pytest.fixture(scope = "session")
def store_raw():
    return SignatureStore.from_parquet(path=config.CRC_SIG, pos_label=None)

@pytest.fixture
def cohort(store):
    m = store.meta
    d = m[m.label == 1].study.value_counts().idxmax()
    n = m[(m.label == 0) & (m.study != d)].study.value_counts().idxmax()
    pids = m.index[((m.study == d) & (m.label == 1)) | ((m.study == n) & (m.label == 0))].tolist()
    return pids, {d, n}


def _studies(store, refs):
    return set(store.meta.loc[refs, "study"])

def _labels(store, refs):
    return store.meta.loc[refs, "label"].tolist()


def test_store_columns(store):
    assert {"label", "label_raw", "study", "assay", "n_cells"} <= set(store.meta.columns)
    assert set(store.meta.label.unique()) == {0, 1}
    assert set(store.meta.label_raw.unique()) == {"tumor", "normal"}

def test_store_sigs_order(store):
    pids = store.meta.index[:5].tolist()
    assert store.sigs(pids).shape == (5, store.sig.shape[1])


def test_cow_leak_free(store, cohort):
    pids, cstud = cohort
    refs = diverse_refs(store, pids, K_dis = 4, K_nor = 4)
    assert not (_studies(store, refs) & cstud)

def test_cow_balance(store, cohort):
    pids, _ = cohort
    labs = _labels(store, diverse_refs(store, pids, K_dis = 4, K_nor = 4))
    assert labs.count(1) == 4 and labs.count(0) == 4

def test_cow_asymmetry(store, cohort):
    pids, _ = cohort
    refs = diverse_refs(store, pids, K_dis = 5, K_nor = 0)
    assert set(_labels(store, refs)) == {1}

def test_cow_no_reuse(store, cohort):
    pids, _ = cohort
    first = diverse_refs(store, pids, K_dis = 3, K_nor = 3)
    second = diverse_refs(store, pids, K_dis = 3, K_nor = 3, used = first)
    assert not (set(first) & set(second))

def test_cow_diversity(store, cohort):
    pids, _ = cohort
    refs = diverse_refs(store, pids, K_dis = 4, K_nor = 4)
    assert len(_studies(store, refs)) > 1


def test_pig_leak_free(store, cohort):
    pids, cstud = cohort
    refs = matched_refs(store, pids, K_case = 4, K_control = 4)
    assert not (_studies(store, refs) & cstud)

def test_pig_cross_labels(store, cohort):
    pids, _ = cohort
    labs = _labels(store, matched_refs(store, pids, K_case = 4, K_control = 4))
    assert labs.count(1) == 4 and labs.count(0) == 4

def test_pig_coverage(store, cohort):
    pids, _ = cohort
    refs = matched_refs(store, pids, K_case = 4, K_control = 4)
    assert len(_studies(store, refs)) > 1

def test_pig_no_reuse(store, cohort):
    pids, _ = cohort
    first = matched_refs(store, pids, K_case = 3, K_control = 3)
    second = matched_refs(store, pids, K_case = 3, K_control = 3, used = first)
    assert not (set(first) & set(second))

def test_pig_underfill(store):
    m = store.meta
    d = m[m.label == 1].study.value_counts().idxmax()
    n = m[(m.label == 0) & (m.study != d)].study.value_counts().idxmax()
    few = m.index[(m.study == d) & (m.label == 1)].tolist()[:3]
    controls = m.index[(m.study == n) & (m.label == 0)].tolist()
    refs = matched_refs(store, few + controls, K_case = 0, K_control = 8)
    assert len(refs) == 8
    assert set(_labels(store, refs)) == {0}
    assert len(_studies(store, refs)) > 1


def test_cow_paired_leak_free(store, cohort):
    pids, cstud = cohort
    refs = diverse_refs_paired(store, pids, K = 3)
    assert refs and not (_studies(store, refs) & cstud)


def test_representatives(store):
    pids = store.meta.index[store.meta.study == store.meta.study.value_counts().idxmax()].tolist()
    reps = representatives(store, pids)
    assert reps and all(p in pids for p in reps)


def test_multiclass_mode(store_raw, cohort):
    pids, _ = cohort
    refs = diverse_refs(store_raw, pids, K_dis = 3, K_nor = 3, case = "tumor", control = "normal")
    assert set(store_raw.meta.loc[refs, "label"]) == {"tumor", "normal"}

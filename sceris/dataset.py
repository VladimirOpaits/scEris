from pathlib import Path
from .store import SignatureStore


def load_cohort(spec, source=None):
    kind = spec.get("kind")
    if kind == "parquet":
        return SignatureStore.from_path(spec["path"])
    if kind == "cells":
        if source is None:
            raise ValueError("raw-cell cohort needs a source to project into")
        return source.embed(spec["path"])
    raise ValueError(f"unknown dataset kind: {kind}")


def spec_for(path):
    return {"kind": "parquet", "path": str(Path(path).resolve())}

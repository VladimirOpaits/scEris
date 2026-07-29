from pathlib import Path
from .store import SignatureStore


class Source:
    def __init__(self, spec):
        self.spec = spec

    def load(self):
        kind = self.spec.get("kind")
        if kind == "local":
            return SignatureStore.from_path(self.spec["path"])
        if kind == "census":
            raise NotImplementedError("census source not wired yet")
        raise ValueError(f"unknown source kind: {kind}")

    def embed(self, cells):
        raise NotImplementedError("cohort projection into source space not wired yet")

    @classmethod
    def local(cls, path):
        return cls({"kind": "local", "path": str(Path(path).resolve())})

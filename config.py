import os
from pathlib import Path

ROOT = Path(__file__).resolve().parent

def _load(name=".env"):
    env = {}
    p = ROOT / name
    if p.exists():
        for line in p.read_text().splitlines():
            line = line.strip()
            if line and not line.startswith("#") and "=" in line:
                k, v = line.split("=", 1)
                env[k.strip()] = v.strip()
    return env

_ENV = _load()

def get(key, default=None):
    return os.environ.get(key, _ENV.get(key, default))   # реальный os.environ перебивает .env

def path(key, default=None):
    v = get(key, default)
    if v is None:
        return None
    return v if os.path.isabs(v) else str(ROOT / v)       # относительные пути -> от корня репо

DATA_DIR   = path("DATA_DIR", "data")
CRC_SIG    = path("CRC_SIG")
CRC_SIG_LF = path("CRC_SIG_LF")
BLOOD_SIG  = path("BLOOD_SIG")
BRAIN_SIG  = path("BRAIN_SIG")

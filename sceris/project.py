import json
from pathlib import Path

STATE = Path(".sceris.json")


def load():
    return json.loads(STATE.read_text()) if STATE.exists() else {}


def save(state):
    STATE.write_text(json.dumps(state, indent=2))


def get(key):
    return load().get(key)


def set(key, value):
    state = load()
    state[key] = value
    save(state)
    return state

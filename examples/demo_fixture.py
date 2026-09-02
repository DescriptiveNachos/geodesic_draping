from __future__ import annotations

import json
from pathlib import Path

import numpy as np


def load_demo_part():
    fixture_dir = Path(__file__).resolve().parents[1] / "test_data" / "fixtures" / "demo_part"
    mesh = json.loads((fixture_dir / "mesh.json").read_text())
    inputs = json.loads((fixture_dir / "inputs.json").read_text())
    return (
        np.asarray(mesh["vertices"], dtype=float),
        np.asarray(mesh["faces"], dtype=np.int64),
        np.asarray(inputs["seed_xy"], dtype=float),
        float(inputs["angle_degrees"]),
    )

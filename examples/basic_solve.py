from __future__ import annotations

import json
from pathlib import Path

import numpy as np

import geodesic_draping as gd


def load_demo_part():
    fixture_dir = Path(__file__).resolve().parents[1] / "test_data" / "fixtures" / "demo_part"
    mesh = json.loads((fixture_dir / "mesh.json").read_text())
    inputs = json.loads((fixture_dir / "inputs.json").read_text())
    vertices = np.asarray(mesh["vertices"], dtype=float)
    faces = np.asarray(mesh["faces"], dtype=np.int64)
    seed_xy = np.asarray(inputs["seed_xy"], dtype=float)
    fabric_angle = float(inputs["angle_degrees"])
    return vertices, faces, seed_xy, fabric_angle


def main() -> None:
    vertices, faces, seed_xy, fabric_angle = load_demo_part()

    result = gd.solve_drape(
        vertices,
        faces,
        seed_xy,
        fabric_angle,
        mode="complete",
        sample_vertex_shear=True,
    )

    print(f"domain: {result.domain}")
    print(f"mode: {result.mode}")
    print(f"vertices: {result.vertices.shape}")
    print(f"faces: {result.faces.shape}")
    print(f"distances: {result.distances.shape if result.distances is not None else None}")
    print(f"vertex shear: {result.vertex_shear.shape if result.vertex_shear is not None else None}")


if __name__ == "__main__":
    main()

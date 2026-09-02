from __future__ import annotations

from pathlib import Path
import sys

import numpy as np
import pytest

import geodesic_draping
from geodesic_draping import GeoDrapeSolver, solve_drape

sys.path.append(str(Path(__file__).resolve().parents[2] / "examples"))
from demo_fixture import load_demo_part


def test_public_api_surface():
    assert geodesic_draping.__version__ == "0.1.0"
    assert geodesic_draping.__all__ == [
        "DrapeResult",
        "GeoDrapeSolver",
        "__version__",
        "solve_drape",
    ]


def test_one_shot_fast_subdivision_shapes():
    vertices, faces, seed_xy, fabric_angle = load_demo_part()

    result = solve_drape(vertices, faces, seed_xy, fabric_angle, mode="fast", retrieval="subdivision")

    assert result.domain == "subdivision"
    assert result.mode == "fast"
    assert result.vertices.shape[1] == 3
    assert result.faces.shape[1] == 3
    assert result.direction_fields.shape == (2, result.faces.shape[0], 3)
    assert result.face_shear.shape == (result.faces.shape[0],)
    assert result.distances is None
    assert result.vertex_shear is None
    assert len(result.generators) == 2
    assert all(len(family) == 2 for family in result.generators)
    assert np.isfinite(result.face_shear).all()


def test_persistent_complete_and_subdivision_retrieve():
    vertices, faces, seed_xy, fabric_angle = load_demo_part()
    solver = GeoDrapeSolver(vertices, faces, intrinsic_backend="integer")

    result = solver.solve(seed_xy, fabric_angle, mode="complete", sample_vertex_shear=True)

    assert result.domain == "extrinsic"
    assert result.mode == "complete"
    assert result.direction_fields is None
    assert result.distances.shape == (2, vertices.shape[0])
    assert result.face_shear is None
    assert result.vertex_shear.shape == (vertices.shape[0],)

    subdivision = solver.retrieve(retrieval="subdivision", sample_vertex_shear=True)

    assert subdivision.domain == "subdivision"
    assert subdivision.mode == "complete"
    assert subdivision.vertices.shape[1] == 3
    assert subdivision.faces.shape[1] == 3
    assert subdivision.direction_fields.shape[0] == 2
    assert subdivision.distances.shape[0] == 2
    assert subdivision.face_shear.shape == (subdivision.faces.shape[0],)


def test_option_validation():
    vertices, faces, seed_xy, fabric_angle = load_demo_part()
    solver = GeoDrapeSolver(vertices, faces)

    with pytest.raises(ValueError):
        solver.solve(seed_xy, fabric_angle, mode="surprise")

    with pytest.raises(NotImplementedError):
        solver.solve(seed_xy, fabric_angle, retrieval="intrinsic")

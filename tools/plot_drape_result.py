from __future__ import annotations
import json
from pathlib import Path
import numpy as np
import geodesic_draping as gd

def plot_result(result, *, show: bool = True):
    try:
        import polyscope as ps
    except ImportError as exc:
        raise ImportError(
            "Python Polyscope is required for tools/plot_drape_result.py. "
            "Install with: python -m pip install polyscope"
        ) from exc

    ps.init()
    mesh = ps.register_surface_mesh("drape result", result.vertices, result.faces)

    if result.face_shear is not None:
        mesh.add_scalar_quantity("face shear", result.face_shear, defined_on="faces", enabled=True)
    if result.vertex_shear is not None:
        mesh.add_scalar_quantity("vertex shear", result.vertex_shear, defined_on="vertices")
    if result.distances is not None:
        mesh.add_scalar_quantity("distance 0", result.distances[0], defined_on="vertices")
        mesh.add_scalar_quantity("distance 1", result.distances[1], defined_on="vertices")
    if result.direction_fields is not None:
        mesh.add_vector_quantity("direction field 0", result.direction_fields[0], defined_on="faces")
        mesh.add_vector_quantity("direction field 1", result.direction_fields[1], defined_on="faces")

    for family_index, family in enumerate(result.generators):
        for direction_index, points in enumerate(family):
            points = np.asarray(points, dtype=float)
            edges = np.array([[i,i+1] for i in range(len(points)-1)])
            if points.shape[0] >= 2:
                ps.register_curve_network(
                    f"generator {family_index}.{direction_index}",
                    points,
                    edges
                )

    if show:
        ps.show()
    return mesh


def _load_fixture(fixture_name):
    fixture_dir = Path(__file__).resolve().parents[1] / "test_data" / "fixtures" / fixture_name
    mesh = json.loads((fixture_dir / "mesh.json").read_text())
    inputs = json.loads((fixture_dir / "inputs.json").read_text())
    return mesh,inputs


if __name__ == "__main__":
    mesh, inputs = _load_fixture('demo_part')

    vertices = np.asarray(mesh["vertices"], dtype=float)
    faces = np.asarray(mesh["faces"], dtype=np.int64)
    seed_xy = np.asarray(inputs["seed_xy"], dtype=float)
    fabric_angle = float(inputs["angle_degrees"])
    
    solver = gd.GeoDrapeSolver(vertices,faces,intrinsic_backend='signpost',refinement='flip')
    result = solver.solve(seed_xy,fabric_angle,mode='complete',retrieval='subdivision',sample_vertex_shear=True) 

    plot_result(result)

from __future__ import annotations

import numpy as np

from basic_solve import load_demo_part

import geodesic_draping as gd


def register_generators(ps, generators):
    for family_index, family in enumerate(generators):
        for direction_index, points in enumerate(family):
            points = np.asarray(points, dtype=float)
            if points.shape[0] < 2:
                continue
            edges = np.array([[i, i + 1] for i in range(points.shape[0] - 1)], dtype=np.int64)
            ps.register_curve_network(f"generator {family_index}.{direction_index}", points, edges)


def main() -> None:
    vertices, faces, seed_xy, fabric_angle = load_demo_part()

    solver = gd.GeoDrapeSolver(vertices, faces)
    solver.solve(seed_xy, fabric_angle, mode="complete")
    result = solver.retrieve(retrieval="subdivision")

    try:
        import polyscope as ps
    except ImportError as exc:
        raise ImportError("Install Python Polyscope with: python -m pip install polyscope") from exc

    ps.init()
    mesh = ps.register_surface_mesh("geodesic draping subdivision", result.vertices, result.faces)
    mesh.add_scalar_quantity("face shear", result.face_shear, defined_on="faces", enabled=True)
    mesh.add_vector_quantity("direction field 0", result.direction_fields[0], defined_on="faces")
    mesh.add_vector_quantity("direction field 1", result.direction_fields[1], defined_on="faces")
    register_generators(ps, result.generators)
    ps.show()


if __name__ == "__main__":
    main()

from __future__ import annotations

from basic_solve import load_demo_part

import geodesic_draping as gd


def main() -> None:
    vertices, faces, seed_xy, fabric_angle = load_demo_part()

    solver = gd.GeoDrapeSolver(vertices, faces, intrinsic_backend="integer")

    fast = solver.solve(seed_xy, fabric_angle, mode="fast", sample_vertex_shear=True)
    complete = solver.solve(seed_xy, fabric_angle, mode="complete", sample_vertex_shear=True)
    subdivision = solver.retrieve(retrieval="subdivision")

    print(f"fast vertex shear: {fast.vertex_shear.shape}")
    print(f"complete distances: {complete.distances.shape}")
    print(f"subdivision face shear: {subdivision.face_shear.shape}")


if __name__ == "__main__":
    main()

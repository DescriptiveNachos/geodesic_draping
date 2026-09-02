from __future__ import annotations

import geodesic_draping as gd

from demo_fixture import load_demo_part


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

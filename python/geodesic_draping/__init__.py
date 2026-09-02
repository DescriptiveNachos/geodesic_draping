from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np

from . import _core

__version__ = "0.1.0"


@dataclass
class DrapeResult:
    """Result returned by geodesic draping solves.

    Fields which are not produced by the selected mode or retrieval domain are
    set to None.
    """

    vertices: np.ndarray
    faces: np.ndarray
    domain: str
    mode: str
    generators: list[list[np.ndarray]]
    direction_fields: np.ndarray | None = None
    distances: np.ndarray | None = None
    face_shear: np.ndarray | None = None
    vertex_shear: np.ndarray | None = None

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> "DrapeResult":
        return cls(**data)

    def to_dict(self) -> dict[str, Any]:
        return vars(self).copy()


class GeoDrapeSolver:
    """Persistent solver for repeated geodesic draping solves on one mesh."""

    def __init__(
        self,
        vertices: np.ndarray,
        faces: np.ndarray,
        *,
        intrinsic_backend: str = "signpost",
        refinement: str = "none",
        angle_threshold: float | None = None,
        circumradius_threshold: float | None = None,
        max_insertions: int | None = None,
        preserve_source_normals: bool = False,
        level_set_constraint: str = "none",
        soft_level_set_weight: float = -1.0,
        diffusion_time_coefficient: float = 1.0,
    ) -> None:
        """Create a solver from triangle mesh arrays.

        Parameters
        ----------
        vertices:
            Float array with shape ``(V, 3)``.
        faces:
            Integer array with shape ``(F, 3)``.
        intrinsic_backend:
            ``"signpost"`` or ``"integer"``.
        refinement:
            ``"none"``, ``"flip"``, or ``"refine"``.
        """
        self._solver = _core.GeoDrapeSolver(
            vertices,
            faces,
            intrinsic_backend=intrinsic_backend,
            refinement=refinement,
            angle_threshold=angle_threshold,
            circumradius_threshold=circumradius_threshold,
            max_insertions=max_insertions,
            preserve_source_normals=preserve_source_normals,
            level_set_constraint=level_set_constraint,
            soft_level_set_weight=soft_level_set_weight,
            diffusion_time_coefficient=diffusion_time_coefficient,
        )

    def solve(
        self,
        seed_xy: np.ndarray,
        fabric_angle: float,
        *,
        mode: str = "complete",
        fiber_angle: float = 90.0,
        retrieval: str = "extrinsic",
        sample_vertex_shear: bool = False,
        trace_length: float | None = None,
        max_trace_iterations: int | None = None,
    ) -> DrapeResult:
        """Run one solve and return data in the requested retrieval domain.

        Public angles are in degrees. Python v1 supports
        ``retrieval="extrinsic"`` and ``retrieval="subdivision"``.
        """
        return DrapeResult.from_dict(
            self._solver.solve(
                seed_xy,
                fabric_angle,
                mode=mode,
                fiber_angle=fiber_angle,
                retrieval=retrieval,
                sample_vertex_shear=sample_vertex_shear,
                trace_length=trace_length,
                max_trace_iterations=max_trace_iterations,
            )
        )

    def retrieve(
        self,
        *,
        retrieval: str = "extrinsic",
        sample_vertex_shear: bool = False,
    ) -> DrapeResult:
        """Retrieve the last solve in another domain."""
        return DrapeResult.from_dict(
            self._solver.retrieve(
                retrieval=retrieval,
                sample_vertex_shear=sample_vertex_shear,
            )
        )


def solve_drape(
    vertices: np.ndarray,
    faces: np.ndarray,
    seed_xy: np.ndarray,
    fabric_angle: float,
    *,
    mode: str = "complete",
    fiber_angle: float = 90.0,
    retrieval: str = "extrinsic",
    sample_vertex_shear: bool = False,
) -> DrapeResult:
    """Convenience one-shot solve.

    Use ``GeoDrapeSolver`` directly when solving repeatedly on the same mesh.
    """
    solver = GeoDrapeSolver(vertices, faces)
    return solver.solve(
        seed_xy,
        fabric_angle,
        mode=mode,
        fiber_angle=fiber_angle,
        retrieval=retrieval,
        sample_vertex_shear=sample_vertex_shear,
    )


__all__ = ["DrapeResult", "GeoDrapeSolver", "__version__", "solve_drape"]

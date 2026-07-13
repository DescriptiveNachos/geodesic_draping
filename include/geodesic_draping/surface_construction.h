#pragma once

#include "geodesic_draping/mesh.h"

#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/vertex_position_geometry.h"

#include <memory>

namespace geodesic_draping {

struct GeometryCentralSurface {
  std::unique_ptr<geometrycentral::surface::ManifoldSurfaceMesh> mesh;
  std::unique_ptr<geometrycentral::surface::VertexPositionGeometry> geometry;
};

GeometryCentralSurface makeGeometryCentralSurface(const SurfaceMeshData& meshData);

} // namespace geodesic_draping

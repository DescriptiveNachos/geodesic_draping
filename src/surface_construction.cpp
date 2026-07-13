#include "geodesic_draping/surface_construction.h"

#include <Eigen/Dense>

#include <limits>
#include <stdexcept>

namespace geodesic_draping {

GeometryCentralSurface makeGeometryCentralSurface(const SurfaceMeshData& meshData) {
  if (meshData.vertices.empty()) {
    throw std::runtime_error("cannot construct geometry-central surface from an empty vertex list");
  }
  if (meshData.faces.empty()) {
    throw std::runtime_error("cannot construct geometry-central surface from an empty face list");
  }
  if (meshData.vertices.size() > static_cast<size_t>(std::numeric_limits<int>::max())) {
    throw std::runtime_error("geometry-central Eigen face matrix adapter currently requires int-sized indices");
  }

  Eigen::MatrixXi faces(static_cast<Eigen::Index>(meshData.faces.size()), 3);
  for (size_t i = 0; i < meshData.faces.size(); ++i) {
    for (size_t j = 0; j < 3; ++j) {
      if (meshData.faces[i][j] >= meshData.vertices.size()) {
        throw std::runtime_error("face index out of range while constructing geometry-central surface");
      }
      faces(static_cast<Eigen::Index>(i), static_cast<Eigen::Index>(j)) =
          static_cast<int>(meshData.faces[i][j]);
    }
  }

  Eigen::MatrixXd vertices(static_cast<Eigen::Index>(meshData.vertices.size()), 3);
  for (size_t i = 0; i < meshData.vertices.size(); ++i) {
    vertices(static_cast<Eigen::Index>(i), 0) = meshData.vertices[i].x();
    vertices(static_cast<Eigen::Index>(i), 1) = meshData.vertices[i].y();
    vertices(static_cast<Eigen::Index>(i), 2) = meshData.vertices[i].z();
  }

  GeometryCentralSurface surface;
  surface.mesh = std::make_unique<geometrycentral::surface::ManifoldSurfaceMesh>(faces);
  surface.geometry =
      std::make_unique<geometrycentral::surface::VertexPositionGeometry>(*surface.mesh, vertices);
  return surface;
}

} // namespace geodesic_draping

#pragma once

#include <Eigen/Dense>

#include <array>
#include <vector>

namespace geodesic_draping {

using Vec2 = Eigen::Vector2d;
using Vec3 = Eigen::Vector3d;
using Face = std::array<size_t, 3>;

struct SurfaceMeshData {
  std::vector<Vec3> vertices;
  std::vector<Face> faces;
};

} // namespace geodesic_draping

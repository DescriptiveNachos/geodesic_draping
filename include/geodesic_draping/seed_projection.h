#pragma once

#include "geodesic_draping/mesh.h"

#include <array>
#include <optional>

namespace geodesic_draping {

struct BarycentricPoint {
  size_t faceIndex = 0;
  Vec3 barycentric = Vec3::Zero();
};

struct SeedProjection {
  Vec3 cartesian = Vec3::Zero();
  BarycentricPoint surfacePoint;
};

std::optional<Vec3> pointInTriangleXY(const Vec2& point, const std::array<Vec3, 3>& triangle);

std::optional<SeedProjection> projectPointXYToMesh(const SurfaceMeshData& mesh,
                                                   const Vec2& point,
                                                   size_t candidateCount = 12);

std::array<Vec3, 4> generateOrthogonalDirections(double angleDegrees,
                                                  const Vec2& seedDirection = Vec2(1.0, 0.0));

} // namespace geodesic_draping

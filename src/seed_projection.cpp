#include "geodesic_draping/seed_projection.h"

#include <algorithm>
#include <cmath>
#include <numeric>

namespace geodesic_draping {
namespace {

std::array<Vec3, 3> triangleForFace(const SurfaceMeshData& mesh, const Face& face) {
  return {mesh.vertices[face[0]], mesh.vertices[face[1]], mesh.vertices[face[2]]};
}

Vec2 centroidXY(const std::array<Vec3, 3>& triangle) {
  return (triangle[0].head<2>() + triangle[1].head<2>() + triangle[2].head<2>()) / 3.0;
}

} // namespace

std::optional<Vec3> pointInTriangleXY(const Vec2& point, const std::array<Vec3, 3>& triangle) {
  const double x = point.x();
  const double y = point.y();
  const double x0 = triangle[0].x();
  const double y0 = triangle[0].y();
  const double x1 = triangle[1].x();
  const double y1 = triangle[1].y();
  const double x2 = triangle[2].x();
  const double y2 = triangle[2].y();

  const double denom = (y1 - y2) * (x0 - x2) + (x2 - x1) * (y0 - y2);
  if (denom == 0.0) {
    return std::nullopt;
  }

  const double u = ((y1 - y2) * (x - x2) + (x2 - x1) * (y - y2)) / denom;
  const double v = ((y2 - y0) * (x - x2) + (x0 - x2) * (y - y2)) / denom;
  const double w = 1.0 - u - v;

  if (u >= 0.0 && v >= 0.0 && w >= 0.0) {
    return Vec3(u, v, w);
  }
  return std::nullopt;
}

std::optional<SeedProjection> projectPointXYToMesh(const SurfaceMeshData& mesh,
                                                   const Vec2& point,
                                                   size_t candidateCount) {
  if (mesh.faces.empty() || mesh.vertices.empty() || candidateCount == 0) {
    return std::nullopt;
  }

  std::vector<size_t> order(mesh.faces.size());
  std::iota(order.begin(), order.end(), size_t{0});

  std::vector<double> centroidDistances(mesh.faces.size(), 0.0);
  for (size_t i = 0; i < mesh.faces.size(); ++i) {
    centroidDistances[i] = (centroidXY(triangleForFace(mesh, mesh.faces[i])) - point).squaredNorm();
  }

  std::sort(order.begin(), order.end(), [&](size_t a, size_t b) {
    if (centroidDistances[a] == centroidDistances[b]) {
      return a < b;
    }
    return centroidDistances[a] < centroidDistances[b];
  });

  const size_t nCandidates = std::min(candidateCount, order.size());
  for (size_t i = 0; i < nCandidates; ++i) {
    const size_t faceIndex = order[i];
    const auto triangle = triangleForFace(mesh, mesh.faces[faceIndex]);
    std::optional<Vec3> barycentric = pointInTriangleXY(point, triangle);
    if (!barycentric) {
      continue;
    }

    SeedProjection projection;
    projection.surfacePoint.faceIndex = faceIndex;
    projection.surfacePoint.barycentric = *barycentric;
    projection.cartesian = (*barycentric)(0) * triangle[0] + (*barycentric)(1) * triangle[1] +
                           (*barycentric)(2) * triangle[2];
    return projection;
  }

  return std::nullopt;
}

std::array<Vec3, 4> generateOrthogonalDirections(double angleDegrees, const Vec2& seedDirection) {
  constexpr double pi = 3.141592653589793238462643383279502884;
  const double angleRadians = angleDegrees * pi / 180.0;
  const double c = std::cos(angleRadians);
  const double s = std::sin(angleRadians);

  const Vec2 vecB(c * seedDirection.x() - s * seedDirection.y(),
                  s * seedDirection.x() + c * seedDirection.y());
  const Vec2 vecC(-vecB.y(), vecB.x());

  const Vec3 b(vecB.x(), vecB.y(), 0.0);
  const Vec3 cVec(vecC.x(), vecC.y(), 0.0);
  return {b, -b, cVec, -cVec};
}

} // namespace geodesic_draping

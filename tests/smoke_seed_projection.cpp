#include "geodesic_draping/seed_projection.h"

#include <cassert>
#include <cmath>

namespace {

bool near(double a, double b, double tolerance = 1e-12) {
  return std::abs(a - b) <= tolerance;
}

void requireBarycentricNear(const geodesic_draping::Vec3& actual,
                            const geodesic_draping::Vec3& expected) {
  assert(near(actual.x(), expected.x()));
  assert(near(actual.y(), expected.y()));
  assert(near(actual.z(), expected.z()));
  assert(near(actual.sum(), 1.0));
}

} // namespace

int main() {
  const std::array<geodesic_draping::Vec3, 3> triangle{
      geodesic_draping::Vec3(0.0, 0.0, 0.0),
      geodesic_draping::Vec3(1.0, 0.0, 0.0),
      geodesic_draping::Vec3(0.0, 1.0, 0.0),
  };

  const auto exactEdge =
      geodesic_draping::pointInTriangleXY(geodesic_draping::Vec2(0.5, 0.0), triangle);
  assert(exactEdge);
  requireBarycentricNear(*exactEdge, geodesic_draping::Vec3(0.5, 0.5, 0.0));

  const auto tinyOutsideEdge =
      geodesic_draping::pointInTriangleXY(geodesic_draping::Vec2(0.5, -1e-14), triangle);
  assert(tinyOutsideEdge);
  requireBarycentricNear(*tinyOutsideEdge, geodesic_draping::Vec3(0.5, 0.5, 0.0));

  const auto clearlyOutsideEdge =
      geodesic_draping::pointInTriangleXY(geodesic_draping::Vec2(0.5, -1e-6), triangle);
  assert(!clearlyOutsideEdge);

  const auto exactVertex =
      geodesic_draping::pointInTriangleXY(geodesic_draping::Vec2(1.0, 0.0), triangle);
  assert(exactVertex);
  requireBarycentricNear(*exactVertex, geodesic_draping::Vec3(0.0, 1.0, 0.0));

  return 0;
}

#include "geodesic_draping/seed_projection.h"
#include "fixture_io.h"

#include <cassert>
#include <cmath>
#include <filesystem>
#include <iostream>
#include <string>
#include <vector>

namespace {

using geodesic_draping::Face;
using geodesic_draping::SurfaceMeshData;
using geodesic_draping::Vec2;
using geodesic_draping::Vec3;

bool near(double a, double b, double tolerance = 1e-10) {
  return std::abs(a - b) <= tolerance;
}

void requireNear(const Vec3& a, const Vec3& b, double tolerance, const std::string& label) {
  if ((a - b).cwiseAbs().maxCoeff() > tolerance) {
    std::cerr << label << " mismatch\n"
              << "actual:   " << a.transpose() << "\n"
              << "expected: " << b.transpose() << "\n";
    assert(false);
  }
}

void testFixture(const std::filesystem::path& root, const std::string& name) {
  const std::filesystem::path fixtureDir = root / name;
  const SurfaceMeshData mesh = geodesic_draping::fixture_io::loadMesh(fixtureDir);
  const Vec2 seedXY = geodesic_draping::fixture_io::loadSeedXY(fixtureDir);
  const double angleDegrees = geodesic_draping::fixture_io::loadAngleDegrees(fixtureDir);

  const auto projection = geodesic_draping::projectPointXYToMesh(mesh, seedXY);
  assert(projection.has_value());
  requireNear(projection->cartesian, geodesic_draping::fixture_io::loadGoldenOrigin(fixtureDir), 1e-9,
              name + " origin");
  assert(near(projection->surfacePoint.barycentric.sum(), 1.0, 1e-10));
  assert(projection->surfacePoint.barycentric.minCoeff() >= -1e-12);

  const auto directions = geodesic_draping::generateOrthogonalDirections(angleDegrees);
  const std::vector<Vec3> goldenDirections = geodesic_draping::fixture_io::loadGoldenDirections(fixtureDir);
  assert(goldenDirections.size() == directions.size());
  for (size_t i = 0; i < directions.size(); ++i) {
    requireNear(directions[i], goldenDirections[i], 1e-12, name + " direction " + std::to_string(i));
  }
}

} // namespace

int main() {
  const std::filesystem::path fixtureRoot = GEODESIC_DRAPING_TEST_DATA_DIR;
  testFixture(fixtureRoot, "tiny_planar");
  testFixture(fixtureRoot, "small_curved");
  testFixture(fixtureRoot, "demo_part");
  return 0;
}

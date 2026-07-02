#include "fixture_io.h"
#include "geodesic_draping/diagnostics.h"

#include <array>
#include <cassert>
#include <cmath>
#include <filesystem>
#include <limits>
#include <stdexcept>

namespace {

void requireNear(double actual, double expected, double tolerance) {
  if (std::abs(actual - expected) > tolerance) {
    throw std::runtime_error("diagnostic value outside tolerance");
  }
}

geodesic_draping::SignedHeatSolveOptions fixtureHeatOptions() {
  geodesic_draping::SignedHeatSolveOptions options;
  options.preserveSourceNormals = false;
  options.levelSetConstraint = geometrycentral::LevelSetConstraint::None;
  options.softLevelSetWeight = -1.0;
  return options;
}

void testSimpleVectors() {
  using geodesic_draping::Vec3;
  const auto diagnostics = geodesic_draping::analyzeVectorMagnitudes({
      Vec3{1.0, 0.0, 0.0},
      Vec3{0.0, 2.0, 0.0},
      Vec3{0.0, 0.0, 0.0},
      Vec3{std::numeric_limits<double>::quiet_NaN(), 0.0, 0.0},
  });

  assert(diagnostics.magnitudes.size() == 4);
  assert(diagnostics.absDeviationFromUnit.size() == 4);
  requireNear(diagnostics.magnitudes[0], 1.0, 1e-15);
  requireNear(diagnostics.magnitudes[1], 2.0, 1e-15);
  requireNear(diagnostics.magnitudes[2], 0.0, 1e-15);
  assert(!std::isfinite(diagnostics.magnitudes[3]));
  requireNear(diagnostics.stats.min, 0.0, 1e-15);
  requireNear(diagnostics.stats.max, 2.0, 1e-15);
  requireNear(diagnostics.stats.mean, 1.0, 1e-15);
  requireNear(diagnostics.stats.maxAbsDeviationFromUnit, 1.0, 1e-15);
  requireNear(diagnostics.stats.meanAbsDeviationFromUnit, 2.0 / 3.0, 1e-15);
  assert(diagnostics.stats.finiteCount == 3);
  assert(diagnostics.stats.nonFiniteCount == 1);
  assert(diagnostics.stats.nearZeroCount == 1);
}

void testFixtureSmoke(const std::filesystem::path& fixtureRoot) {
  const std::filesystem::path fixtureDir = fixtureRoot / "tiny_planar";
  const auto mesh = geodesic_draping::fixture_io::loadMesh(fixtureDir);
  const auto seedXY = geodesic_draping::fixture_io::loadSeedXY(fixtureDir);
  const double angleDegrees = geodesic_draping::fixture_io::loadAngleDegrees(fixtureDir);
  const auto options = fixtureHeatOptions();

  const auto complete = geodesic_draping::solveCompleteDrape(mesh, seedXY, angleDegrees, options);
  const auto fast = geodesic_draping::solveFastDrape(mesh, seedXY, angleDegrees, options);
  const auto completeDiagnostics = geodesic_draping::analyzeCompleteGradientMagnitudes(complete);
  auto surface = geodesic_draping::makeGeometryCentralSurface(mesh);
  const std::array<geodesic_draping::VectorMagnitudeDiagnostics, 2> fastDiagnostics{
      geodesic_draping::analyzeVectorMagnitudes(
          geodesic_draping::faceDirectionsToExtrinsicVectors(surface, fast.faceDirections[0])),
      geodesic_draping::analyzeVectorMagnitudes(
          geodesic_draping::faceDirectionsToExtrinsicVectors(surface, fast.faceDirections[1])),
  };

  for (const auto& diagnostic : completeDiagnostics) {
    assert(diagnostic.magnitudes.size() == mesh.vertices.size());
    assert(diagnostic.stats.finiteCount == mesh.vertices.size());
    assert(std::isfinite(diagnostic.stats.maxAbsDeviationFromUnit));
  }
  for (const auto& diagnostic : fastDiagnostics) {
    assert(diagnostic.magnitudes.size() == mesh.faces.size());
    assert(diagnostic.stats.finiteCount == mesh.faces.size());
    assert(std::isfinite(diagnostic.stats.maxAbsDeviationFromUnit));
  }
}

} // namespace

int main() {
  testSimpleVectors();
  testFixtureSmoke(GEODESIC_DRAPING_TEST_DATA_DIR);
  return 0;
}

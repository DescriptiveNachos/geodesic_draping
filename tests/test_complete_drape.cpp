#include "fixture_io.h"
#include "geodesic_draping/geodrape.h"

#include <cassert>
#include <cmath>
#include <filesystem>
#include <iostream>
#include <string>

namespace {

geodesic_draping::SignedHeatSolveOptions fixtureHeatOptions() {
  geodesic_draping::SignedHeatSolveOptions options;
  options.preserveSourceNormals = false;
  options.levelSetConstraint = geometrycentral::LevelSetConstraint::None;
  options.softLevelSetWeight = -1.0;
  return options;
}

void requireNearArray(const std::vector<double>& actual,
                      const std::vector<double>& golden,
                      double tolerance,
                      const std::string& label) {
  assert(actual.size() == golden.size());
  double maxDiff = 0.0;
  for (size_t i = 0; i < actual.size(); ++i) {
    maxDiff = std::max(maxDiff, std::abs(actual[i] - golden[i]));
  }
  if (maxDiff > tolerance) {
    std::cerr << label << " max diff " << maxDiff << " exceeds tolerance " << tolerance << "\n";
    assert(false);
  }
}

void requireFinite(const geodesic_draping::CompleteDrapeResult& result, size_t nVertices) {
  for (const auto& distance : result.distances) {
    assert(distance.size() == nVertices);
    for (double value : distance) {
      assert(std::isfinite(value));
    }
  }
  for (const auto& gradient : result.gradients) {
    assert(gradient.size() == nVertices);
    for (const auto& value : gradient) {
      assert(std::isfinite(value.x()));
      assert(std::isfinite(value.y()));
      assert(std::isfinite(value.z()));
    }
  }
  assert(result.shearAnglesDegrees.size() == nVertices);
  for (double value : result.shearAnglesDegrees) {
    assert(std::isfinite(value));
  }
}

geodesic_draping::CompleteDrapeResult solveFixture(const std::filesystem::path& root, const std::string& name) {
  const std::filesystem::path fixtureDir = root / name;
  return geodesic_draping::solveCompleteDrape(
      geodesic_draping::fixture_io::loadMesh(fixtureDir),
      geodesic_draping::fixture_io::loadSeedXY(fixtureDir),
      geodesic_draping::fixture_io::loadAngleDegrees(fixtureDir),
      fixtureHeatOptions());
}

void testPersistentSolver(const std::filesystem::path& root, const std::string& name) {
  const std::filesystem::path fixtureDir = root / name;
  geodesic_draping::GeoDrapeSolver solver(
      geodesic_draping::fixture_io::loadMesh(fixtureDir),
      fixtureHeatOptions());
  const auto seedXY = geodesic_draping::fixture_io::loadSeedXY(fixtureDir);
  const double angleDegrees = geodesic_draping::fixture_io::loadAngleDegrees(fixtureDir);

  const auto first = solver.solveComplete(seedXY, angleDegrees);
  const auto second = solver.solveComplete(seedXY, angleDegrees);
  requireNearArray(second.distances[0], first.distances[0], 1e-12, name + " persistent complete dist_0");
  requireNearArray(second.distances[1], first.distances[1], 1e-12, name + " persistent complete dist_1");
  requireNearArray(second.shearAnglesDegrees, first.shearAnglesDegrees, 1e-12, name + " persistent complete shear");
}

void testOneShotComplete(const std::filesystem::path& root, const std::string& name) {
  const std::filesystem::path fixtureDir = root / name;
  const auto mesh = geodesic_draping::fixture_io::loadMesh(fixtureDir);
  const auto seedXY = geodesic_draping::fixture_io::loadSeedXY(fixtureDir);
  const double angleDegrees = geodesic_draping::fixture_io::loadAngleDegrees(fixtureDir);
  geodesic_draping::DrapeSolveOptions solveOptions;
  solveOptions.mode = geodesic_draping::DrapeSolveMode::Complete;

  const auto oneShot = geodesic_draping::solveDrape(mesh, seedXY, angleDegrees, fixtureHeatOptions(), solveOptions);
  const auto complete = geodesic_draping::solveCompleteDrape(mesh, seedXY, angleDegrees, fixtureHeatOptions());
  assert(oneShot.distances);
  assert(oneShot.gradients);
  assert(!oneShot.faceShearAnglesDegrees);
  assert(oneShot.vertexShearAnglesDegrees);
  requireNearArray((*oneShot.distances)[0], complete.distances[0], 1e-12, name + " one-shot complete dist_0");
  requireNearArray((*oneShot.distances)[1], complete.distances[1], 1e-12, name + " one-shot complete dist_1");
  requireNearArray(*oneShot.vertexShearAnglesDegrees,
                   complete.shearAnglesDegrees,
                   1e-12,
                   name + " one-shot complete shear");

  solveOptions.sampleSecondaryShear = true;
  const auto sampled = geodesic_draping::solveDrape(mesh, seedXY, angleDegrees, fixtureHeatOptions(), solveOptions);
  assert(sampled.faceShearAnglesDegrees);
  assert(sampled.faceShearAnglesDegrees->size() == mesh.faces.size());
}

} // namespace

int main() {
  const std::filesystem::path fixtureRoot = GEODESIC_DRAPING_TEST_DATA_DIR;

  const auto tiny = solveFixture(fixtureRoot, "tiny_planar");
  requireFinite(tiny, geodesic_draping::fixture_io::loadMesh(fixtureRoot / "tiny_planar").vertices.size());
  requireNearArray(tiny.distances[0],
                   geodesic_draping::fixture_io::loadGoldenScalarArray(fixtureRoot / "tiny_planar", "dist_0"),
                   1e-8,
                   "tiny complete dist_0");
  requireNearArray(tiny.distances[1],
                   geodesic_draping::fixture_io::loadGoldenScalarArray(fixtureRoot / "tiny_planar", "dist_1"),
                   1e-8,
                   "tiny complete dist_1");
  requireNearArray(tiny.shearAnglesDegrees,
                   geodesic_draping::fixture_io::loadGoldenShearArray(fixtureRoot / "tiny_planar", "complete"),
                   1e-8,
                   "tiny complete shear");

  const auto small = solveFixture(fixtureRoot, "small_curved");
  requireFinite(small, geodesic_draping::fixture_io::loadMesh(fixtureRoot / "small_curved").vertices.size());

  const auto smoothGood = solveFixture(fixtureRoot, "smooth_quality_good");
  requireFinite(smoothGood,
                geodesic_draping::fixture_io::loadMesh(fixtureRoot / "smooth_quality_good").vertices.size());

  const auto smoothPoor = solveFixture(fixtureRoot, "smooth_quality_poor");
  requireFinite(smoothPoor,
                geodesic_draping::fixture_io::loadMesh(fixtureRoot / "smooth_quality_poor").vertices.size());

  const auto demo = solveFixture(fixtureRoot, "demo_part");
  requireFinite(demo, geodesic_draping::fixture_io::loadMesh(fixtureRoot / "demo_part").vertices.size());
  testPersistentSolver(fixtureRoot, "demo_part");
  testOneShotComplete(fixtureRoot, "tiny_planar");

  return 0;
}

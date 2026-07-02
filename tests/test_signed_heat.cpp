#include "fixture_io.h"
#include "geodesic_draping/generator_tracing.h"
#include "geodesic_draping/seed_projection.h"
#include "geodesic_draping/signed_heat.h"
#include "geodesic_draping/custom_signed_heat.h"

#include <array>
#include <cassert>
#include <cmath>
#include <filesystem>
#include <iostream>
#include <string>

namespace {

struct FixtureSolve {
  geodesic_draping::SurfaceMeshData meshData;
  std::array<std::vector<double>, 2> distances;
};

FixtureSolve solveFixture(const std::filesystem::path& root, const std::string& name) {
  const std::filesystem::path fixtureDir = root / name;
  geodesic_draping::SurfaceMeshData meshData = geodesic_draping::fixture_io::loadMesh(fixtureDir);
  auto surface = geodesic_draping::makeGeometryCentralSurface(meshData);

  const geodesic_draping::Vec2 seedXY = geodesic_draping::fixture_io::loadSeedXY(fixtureDir);
  const double angleDegrees = geodesic_draping::fixture_io::loadAngleDegrees(fixtureDir);
  const auto projection = geodesic_draping::projectPointXYToMesh(meshData, seedXY);
  assert(projection.has_value());

  const auto directions = geodesic_draping::generateOrthogonalDirections(angleDegrees);
  const auto traces = geodesic_draping::traceGenerators(surface, projection->surfacePoint, directions);
  const auto sourceCurves = geodesic_draping::pairOppositeGeneratorTraces(traces);

  geodesic_draping::SignedHeatSolveOptions options;
  options.preserveSourceNormals = false;
  options.levelSetConstraint = geometrycentral::LevelSetConstraint::None;
  options.softLevelSetWeight = -1.0;

  geodesic_draping::CustomSignedHeatSolver solver(surface, options.diffusionTimeCoefficient);
  const auto heatSolves = solver.solve(sourceCurves, options, true);
  return {meshData, {heatSolves[0].distance, heatSolves[1].distance}};
}

void requireFiniteDistanceFields(const FixtureSolve& solve, const std::string& name) {
  for (size_t fieldIndex = 0; fieldIndex < solve.distances.size(); ++fieldIndex) {
    const auto& field = solve.distances[fieldIndex];
    assert(field.size() == solve.meshData.vertices.size());
    for (double value : field) {
      if (!std::isfinite(value)) {
        std::cerr << name << " distance field " << fieldIndex << " contains non-finite value\n";
        assert(false);
      }
    }
  }
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

} // namespace

int main() {
  const std::filesystem::path fixtureRoot = GEODESIC_DRAPING_TEST_DATA_DIR;

  const FixtureSolve tiny = solveFixture(fixtureRoot, "tiny_planar");
  requireFiniteDistanceFields(tiny, "tiny_planar");
  requireNearArray(tiny.distances[0],
                   geodesic_draping::fixture_io::loadGoldenScalarArray(fixtureRoot / "tiny_planar", "dist_0"),
                   1e-8,
                   "tiny_planar dist_0");
  requireNearArray(tiny.distances[1],
                   geodesic_draping::fixture_io::loadGoldenScalarArray(fixtureRoot / "tiny_planar", "dist_1"),
                   1e-8,
                   "tiny_planar dist_1");

  requireFiniteDistanceFields(solveFixture(fixtureRoot, "small_curved"), "small_curved");
  requireFiniteDistanceFields(solveFixture(fixtureRoot, "demo_part"), "demo_part");

  return 0;
}

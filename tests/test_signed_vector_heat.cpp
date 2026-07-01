#include "fixture_io.h"
#include "geodesic_draping/geodrape.h"
#include "geodesic_draping/geometrycentral_adapter.h"
#include "geodesic_draping/seed_projection.h"
#include "geodesic_draping/signed_vector_heat.h"

#include <cmath>
#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <string>

namespace {

geodesic_draping::SignedHeatSolveOptions fixtureHeatOptions() {
  geodesic_draping::SignedHeatSolveOptions options;
  options.preserveSourceNormals = false;
  options.levelSetConstraint = geometrycentral::LevelSetConstraint::None;
  options.softLevelSetWeight = -1.0;
  return options;
}

void requireNearVectorArray(const std::vector<geodesic_draping::Vec3>& actual,
                            const std::vector<geodesic_draping::Vec3>& golden,
                            double tolerance,
                            const std::string& label) {
  assert(actual.size() == golden.size());
  double maxDiff = 0.0;
  for (size_t i = 0; i < actual.size(); ++i) {
    maxDiff = std::max(maxDiff, (actual[i] - golden[i]).cwiseAbs().maxCoeff());
  }
  if (maxDiff > tolerance) {
    std::cerr << label << " max diff " << maxDiff << " exceeds tolerance " << tolerance << "\n";
    throw std::runtime_error(label + " exceeded tolerance");
  }
}

void testFixture(const std::filesystem::path& root,
                 const std::string& name,
                 double tolerance) {
  const std::filesystem::path fixtureDir = root / name;
  const auto mesh = geodesic_draping::fixture_io::loadMesh(fixtureDir);
  const auto seedXY = geodesic_draping::fixture_io::loadSeedXY(fixtureDir);
  const double angleDegrees = geodesic_draping::fixture_io::loadAngleDegrees(fixtureDir);

  const auto seed = geodesic_draping::projectPointXYToMesh(mesh, seedXY);
  if (!seed) {
    throw std::runtime_error(name + " seed projection failed");
  }
  const auto directions = geodesic_draping::generateOrthogonalDirections(angleDegrees);
  auto surface = geodesic_draping::makeGeometryCentralSurface(mesh);
  const auto generators = geodesic_draping::traceGenerators(surface, seed->surfacePoint, directions);
  const auto sourceCurves = geodesic_draping::pairOppositeGeneratorTraces(generators);
  const auto options = fixtureHeatOptions();
  const auto heat = geodesic_draping::computeSignedVectorHeats(surface, sourceCurves, options);

  requireNearVectorArray(heat[0].vertexVectorHeat,
                         geodesic_draping::fixture_io::loadGoldenVectorArray(fixtureDir, "fast", "grad_0"),
                         tolerance,
                         name + " fast grad_0");
  requireNearVectorArray(heat[1].vertexVectorHeat,
                         geodesic_draping::fixture_io::loadGoldenVectorArray(fixtureDir, "fast", "grad_1"),
                         tolerance,
                         name + " fast grad_1");

  if (heat[0].diffusion.sourceEdgeVectorHeat.size() != surface.mesh->nEdges() ||
      heat[0].diffusion.diffusedEdgeVectorHeat.size() != surface.mesh->nEdges() ||
      heat[0].normalizedFaceVectorHeat.size() != surface.mesh->nFaces()) {
    throw std::runtime_error(name + " fast heat intermediate sizes do not match mesh");
  }

  const auto fast = geodesic_draping::solveFastDrape(mesh, seedXY, angleDegrees, options);
  requireNearVectorArray(fast.gradients[0],
                         geodesic_draping::fixture_io::loadGoldenVectorArray(fixtureDir, "fast", "grad_0"),
                         tolerance,
                         name + " high-level fast grad_0");
  requireNearVectorArray(fast.gradients[1],
                         geodesic_draping::fixture_io::loadGoldenVectorArray(fixtureDir, "fast", "grad_1"),
                         tolerance,
                         name + " high-level fast grad_1");
}

} // namespace

int main() {
  const std::filesystem::path fixtureRoot = GEODESIC_DRAPING_TEST_DATA_DIR;
  testFixture(fixtureRoot, "tiny_planar", 1e-12);
  testFixture(fixtureRoot, "small_curved", 2e-2);
  testFixture(fixtureRoot, "demo_part", 1e-10);
  return 0;
}

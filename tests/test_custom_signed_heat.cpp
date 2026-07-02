#include "fixture_io.h"
#include "geodesic_draping/geodrape.h"
#include "geodesic_draping/geometrycentral_adapter.h"
#include "geodesic_draping/seed_projection.h"
#include "geodesic_draping/custom_signed_heat.h"

#include <array>
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

void requireNearArray(const std::vector<double>& actual,
                      const std::vector<double>& expected,
                      double tolerance,
                      const std::string& label) {
  if (actual.size() != expected.size()) {
    throw std::runtime_error(label + " size mismatch");
  }
  double maxDiff = 0.0;
  for (size_t i = 0; i < actual.size(); ++i) {
    maxDiff = std::max(maxDiff, std::abs(actual[i] - expected[i]));
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
  const auto heat = geodesic_draping::computeCustomSignedHeat(surface, sourceCurves, options);

  if (heat[0].diffusion.sourceEdgeHeatField.size() != surface.mesh->nEdges() ||
      heat[0].diffusion.diffusedEdgeHeatField.size() != surface.mesh->nEdges() ||
      heat[0].normalizedFaceDirections.size() != surface.mesh->nFaces()) {
    throw std::runtime_error(name + " fast heat intermediate sizes do not match mesh");
  }
  const auto faceShear = geodesic_draping::computeFaceShearAnglesDegrees(
      surface,
      heat[0].normalizedFaceDirections,
      heat[1].normalizedFaceDirections);
  if (faceShear.size() != surface.mesh->nFaces()) {
    throw std::runtime_error(name + " face shear size does not match mesh");
  }
  for (double value : faceShear) {
    if (!std::isfinite(value)) {
      throw std::runtime_error(name + " face shear contains non-finite values");
    }
  }

  const auto vertexShear = geodesic_draping::averageFaceScalarsToVertices(mesh, faceShear);
  if (vertexShear.size() != mesh.vertices.size()) {
    throw std::runtime_error(name + " sampled vertex shear does not match mesh");
  }

  const auto fast = geodesic_draping::solveFastDrape(mesh, seedXY, angleDegrees, options);
  if (fast.faceDirections[0].size() != surface.mesh->nFaces() ||
      fast.faceDirections[1].size() != surface.mesh->nFaces() ||
      fast.faceShearAnglesDegrees.size() != surface.mesh->nFaces() ||
      fast.vertexShearAnglesDegrees.size() != mesh.vertices.size()) {
    throw std::runtime_error(name + " high-level fast output sizes do not match mesh");
  }
  if (!fast.distances[0].empty() || !fast.distances[1].empty()) {
    throw std::runtime_error(name + " default fast result should not return distances");
  }
  requireNearArray(fast.faceShearAnglesDegrees, faceShear, tolerance, name + " high-level face shear");
  requireNearArray(fast.vertexShearAnglesDegrees, vertexShear, tolerance, name + " high-level vertex shear");

  geodesic_draping::DrapeSolveOptions hybridOptions;
  hybridOptions.mode = geodesic_draping::DrapeSolveMode::Hybrid;
  const auto fastWithDistances = geodesic_draping::solveFastDrape(mesh, seedXY, angleDegrees, options, hybridOptions);
  const auto hybrid = geodesic_draping::solveDrape(mesh, seedXY, angleDegrees, options, hybridOptions);
  const auto complete = geodesic_draping::solveCompleteDrape(mesh, seedXY, angleDegrees, options);
  requireNearArray(fastWithDistances.distances[0], complete.distances[0], 1e-12, name + " hybrid dist_0");
  requireNearArray(fastWithDistances.distances[1], complete.distances[1], 1e-12, name + " hybrid dist_1");
  requireNearArray(hybrid.distances[0], complete.distances[0], 1e-12, name + " one-shot hybrid dist_0");
  requireNearArray(hybrid.distances[1], complete.distances[1], 1e-12, name + " one-shot hybrid dist_1");
  requireNearArray(hybrid.faceShearAnglesDegrees,
                   fast.faceShearAnglesDegrees,
                   tolerance,
                   name + " one-shot hybrid preserves face shear");
  requireNearArray(fastWithDistances.faceShearAnglesDegrees,
                   fast.faceShearAnglesDegrees,
                   tolerance,
                   name + " hybrid preserves face shear");
  requireNearArray(fastWithDistances.vertexShearAnglesDegrees,
                   fast.vertexShearAnglesDegrees,
                   tolerance,
                   name + " hybrid preserves vertex shear");

  geodesic_draping::GeoDrapeSolver solver(mesh, options);
  const auto first = solver.solveFast(seedXY, angleDegrees);
  const auto second = solver.solveFast(seedXY, angleDegrees);
  requireNearArray(second.faceShearAnglesDegrees, first.faceShearAnglesDegrees, 1e-12, name + " persistent fast face shear");
  requireNearArray(second.vertexShearAnglesDegrees, first.vertexShearAnglesDegrees, 1e-12, name + " persistent fast vertex shear");
}

void testFastFinite(const std::filesystem::path& root, const std::string& name) {
  const std::filesystem::path fixtureDir = root / name;
  const auto mesh = geodesic_draping::fixture_io::loadMesh(fixtureDir);
  const auto seedXY = geodesic_draping::fixture_io::loadSeedXY(fixtureDir);
  const double angleDegrees = geodesic_draping::fixture_io::loadAngleDegrees(fixtureDir);
  const auto fast = geodesic_draping::solveFastDrape(mesh, seedXY, angleDegrees, fixtureHeatOptions());

  if (fast.faceShearAnglesDegrees.size() != mesh.faces.size() ||
      fast.vertexShearAnglesDegrees.size() != mesh.vertices.size()) {
    throw std::runtime_error(name + " fast output sizes do not match mesh");
  }
  for (const auto& field : fast.faceDirections) {
    if (field.size() != mesh.faces.size()) {
      throw std::runtime_error(name + " fast face direction sizes do not match mesh");
    }
    for (const auto& value : field) {
      if (!std::isfinite(value.x()) || !std::isfinite(value.y()) || !std::isfinite(value.z())) {
        throw std::runtime_error(name + " fast face direction field contains non-finite values");
      }
    }
  }
  for (double value : fast.faceShearAnglesDegrees) {
    if (!std::isfinite(value)) {
      throw std::runtime_error(name + " fast face shear contains non-finite values");
    }
  }
  for (double value : fast.vertexShearAnglesDegrees) {
    if (!std::isfinite(value)) {
      throw std::runtime_error(name + " fast vertex shear contains non-finite values");
    }
  }
}

void testAnalyticFaceShear() {
  geodesic_draping::SurfaceMeshData mesh;
  mesh.vertices = {
      {0.0, 0.0, 0.0},
      {1.0, 0.0, 0.0},
      {0.0, 1.0, 0.0},
  };
  mesh.faces = {{{0, 1, 2}}};
  auto surface = geodesic_draping::makeGeometryCentralSurface(mesh);

  const geodesic_draping::FaceHeatDirectionField xDirection{{-1.0, 1.0, 0.0}};
  const geodesic_draping::FaceHeatDirectionField yDirection{{-1.0, 0.0, 1.0}};

  const auto orthogonalShear =
      geodesic_draping::computeFaceShearAnglesDegrees(surface, xDirection, yDirection);
  if (orthogonalShear.size() != 1 || std::abs(orthogonalShear[0]) > 1e-12) {
    throw std::runtime_error("orthogonal analytic face shear should be zero");
  }

  const auto parallelShear =
      geodesic_draping::computeFaceShearAnglesDegrees(surface, xDirection, xDirection);
  if (parallelShear.size() != 1 || std::abs(parallelShear[0] - 90.0) > 1e-12) {
    throw std::runtime_error("parallel analytic face shear should be ninety degrees");
  }

  const geodesic_draping::FaceHeatDirectionField negativeXDirection{{1.0, -1.0, 0.0}};
  const auto oppositeShear =
      geodesic_draping::computeFaceShearAnglesDegrees(surface, xDirection, negativeXDirection);
  if (oppositeShear.size() != 1 || std::abs(oppositeShear[0] - 90.0) > 1e-12) {
    throw std::runtime_error("opposite analytic face shear should be ninety degrees");
  }

  const auto extrinsic = geodesic_draping::faceDirectionsToExtrinsicVectors(surface, xDirection);
  if (extrinsic.size() != 1 || (extrinsic[0] - geodesic_draping::Vec3{1.0, 0.0, 0.0}).norm() > 1e-12) {
    throw std::runtime_error("analytic face direction should convert to the expected extrinsic vector");
  }
}

} // namespace

int main() {
  const std::filesystem::path fixtureRoot = GEODESIC_DRAPING_TEST_DATA_DIR;
  testAnalyticFaceShear();
  testFixture(fixtureRoot, "tiny_planar", 1e-12);
  testFixture(fixtureRoot, "small_curved", 2e-2);
  testFixture(fixtureRoot, "demo_part", 1e-10);
  testFastFinite(fixtureRoot, "smooth_quality_good");
  testFastFinite(fixtureRoot, "smooth_quality_poor");
  return 0;
}

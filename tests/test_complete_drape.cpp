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

geodesic_draping::DrapeSolveOptions completeOptions(bool sampleVertexShear = true) {
  geodesic_draping::DrapeSolveOptions options;
  options.mode = geodesic_draping::DrapeSolveMode::Complete;
  options.sampleVertexShear = sampleVertexShear;
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

double maxAbsDiff(const geodesic_draping::Vec3& a, const geodesic_draping::Vec3& b) {
  return (a - b).cwiseAbs().maxCoeff();
}

void requireNearVector(const geodesic_draping::Vec3& actual,
                       const geodesic_draping::Vec3& expected,
                       double tolerance,
                       const std::string& label) {
  const double diff = maxAbsDiff(actual, expected);
  if (diff > tolerance) {
    std::cerr << label << " max diff " << diff << " exceeds tolerance " << tolerance << "\n"
              << "actual:   " << actual.transpose() << "\n"
              << "expected: " << expected.transpose() << "\n";
    throw std::runtime_error(label + " exceeded tolerance");
  }
}

void requireFiniteComplete(const geodesic_draping::DrapeResult& result, size_t nVertices) {
  assert(result.mode == geodesic_draping::DrapeSolveMode::Complete);
  assert(result.distances);
  assert(result.vertexShear);
  for (const auto& distance : *result.distances) {
    assert(distance.size() == nVertices);
    for (double value : distance) {
      assert(std::isfinite(value));
    }
  }
  assert(result.vertexShear->size() == nVertices);
  for (double value : *result.vertexShear) {
    assert(std::isfinite(value));
  }
}

void testHighLevelGeneratorGoldenEndpoints(const std::filesystem::path& root, const std::string& name) {
  const std::filesystem::path fixtureDir = root / name;
  const auto mesh = geodesic_draping::fixture_io::loadMesh(fixtureDir);
  const auto seedXY = geodesic_draping::fixture_io::loadSeedXY(fixtureDir);
  const double angleDegrees = geodesic_draping::fixture_io::loadAngleDegrees(fixtureDir);
  const std::vector<geodesic_draping::Vec3> goldenEnds =
      geodesic_draping::fixture_io::loadGoldenGeneratorLastPoints(fixtureDir);

  geodesic_draping::DrapeSolveOptions options;
  options.mode = geodesic_draping::DrapeSolveMode::Fast;
  options.retrieval = geodesic_draping::ResultDomain::Extrinsic;
  options.sampleVertexShear = false;
  const auto result = geodesic_draping::solveDrape(
      mesh,
      seedXY,
      angleDegrees,
      fixtureHeatOptions(),
      options);

  const std::array<const geodesic_draping::DrapeTrace*, 4> traces = {
      &result.traces[0].positive,
      &result.traces[0].negative,
      &result.traces[1].positive,
      &result.traces[1].negative,
  };
  if (goldenEnds.size() != traces.size()) {
    throw std::runtime_error(name + " high-level generator golden endpoint count mismatch");
  }
  for (size_t i = 0; i < traces.size(); ++i) {
    if (traces[i]->extrinsicPoints.empty()) {
      throw std::runtime_error(name + " high-level generator " + std::to_string(i) + " returned no points");
    }
    requireNearVector(
        traces[i]->extrinsicPoints.back(),
        goldenEnds[i],
        1e-5,
        name + " high-level generator " + std::to_string(i) + " end");
  }
}

geodesic_draping::DrapeResult solveFixture(const std::filesystem::path& root, const std::string& name) {
  const std::filesystem::path fixtureDir = root / name;
  return geodesic_draping::solveDrape(
      geodesic_draping::fixture_io::loadMesh(fixtureDir),
      geodesic_draping::fixture_io::loadSeedXY(fixtureDir),
      geodesic_draping::fixture_io::loadAngleDegrees(fixtureDir),
      fixtureHeatOptions(),
      completeOptions());
}

void testPersistentSolver(const std::filesystem::path& root, const std::string& name) {
  const std::filesystem::path fixtureDir = root / name;
  geodesic_draping::GeoDrapeSolver solver(
      geodesic_draping::fixture_io::loadMesh(fixtureDir),
      fixtureHeatOptions());
  const auto seedXY = geodesic_draping::fixture_io::loadSeedXY(fixtureDir);
  const double angleDegrees = geodesic_draping::fixture_io::loadAngleDegrees(fixtureDir);

  const auto options = completeOptions();
  const auto first = solver.solve(seedXY, angleDegrees, options);
  const auto second = solver.solve(seedXY, angleDegrees, options);
  assert(first.distances && second.distances);
  assert(first.vertexShear && second.vertexShear);
  requireNearArray((*second.distances)[0], (*first.distances)[0], 1e-12, name + " persistent complete dist_0");
  requireNearArray((*second.distances)[1], (*first.distances)[1], 1e-12, name + " persistent complete dist_1");
  requireNearArray(*second.vertexShear, *first.vertexShear, 1e-12, name + " persistent complete shear");
}

void testOneShotComplete(const std::filesystem::path& root, const std::string& name) {
  const std::filesystem::path fixtureDir = root / name;
  const auto mesh = geodesic_draping::fixture_io::loadMesh(fixtureDir);
  const auto seedXY = geodesic_draping::fixture_io::loadSeedXY(fixtureDir);
  const double angleDegrees = geodesic_draping::fixture_io::loadAngleDegrees(fixtureDir);
  geodesic_draping::DrapeSolveOptions solveOptions = completeOptions();

  const auto oneShot = geodesic_draping::solveDrape(mesh, seedXY, angleDegrees, fixtureHeatOptions(), solveOptions);
  assert(oneShot.distances);
  assert(!oneShot.faceShear);
  assert(oneShot.vertexShear);
  requireNearArray((*oneShot.distances)[0], (*oneShot.distances)[0], 1e-12, name + " one-shot complete dist_0");
  requireNearArray((*oneShot.distances)[1], (*oneShot.distances)[1], 1e-12, name + " one-shot complete dist_1");
  requireNearArray(*oneShot.vertexShear,
                   *oneShot.vertexShear,
                   1e-12,
                   name + " one-shot complete shear");

  geodesic_draping::DrapeSolveOptions intrinsicOptions = completeOptions(false);
  intrinsicOptions.retrieval = geodesic_draping::ResultDomain::Intrinsic;
  const auto intrinsic = geodesic_draping::solveDrape(mesh, seedXY, angleDegrees, fixtureHeatOptions(), intrinsicOptions);
  assert(intrinsic.faceShear);
  assert(intrinsic.faceShear->size() == mesh.faces.size());
}

void requireFiniteFaceShear(const geodesic_draping::DrapeResult& result, const std::string& label) {
  if (!result.faceShear) {
    throw std::runtime_error(label + " did not return face shear");
  }
  if (result.faceShear->empty()) {
    throw std::runtime_error(label + " returned empty face shear");
  }
  for (double value : *result.faceShear) {
    if (!std::isfinite(value)) {
      throw std::runtime_error(label + " contains non-finite face shear");
    }
  }
  for (const auto& family : result.traces) {
    const bool positiveEmpty = family.positive.extrinsicPoints.empty() && family.positive.intrinsicPoints.empty();
    const bool negativeEmpty = family.negative.extrinsicPoints.empty() && family.negative.intrinsicPoints.empty();
    if (positiveEmpty || negativeEmpty) {
      throw std::runtime_error(label + " returned an empty generator trace");
    }
  }
}

void testRefinementSmoke(const std::filesystem::path& root) {
  const std::filesystem::path fixtureDir = root / "small_curved";
  const auto mesh = geodesic_draping::fixture_io::loadMesh(fixtureDir);
  const auto seedXY = geodesic_draping::fixture_io::loadSeedXY(fixtureDir);
  const double angleDegrees = geodesic_draping::fixture_io::loadAngleDegrees(fixtureDir);

  geodesic_draping::DrapeSolveOptions options;
  options.mode = geodesic_draping::DrapeSolveMode::Fast;
  options.retrieval = geodesic_draping::ResultDomain::Intrinsic;
  options.sampleVertexShear = false;

  geodesic_draping::RefinementOptions flip;
  flip.mode = geodesic_draping::RefinementMode::DelaunayFlip;
  geodesic_draping::GeoDrapeSolver flipSolver(mesh, fixtureHeatOptions(), {}, flip);
  requireFiniteFaceShear(flipSolver.solve(seedXY, angleDegrees, options), "delaunay flip smoke");

  geodesic_draping::RefinementOptions refine;
  refine.mode = geodesic_draping::RefinementMode::DelaunayRefine;
  refine.maxInsertions = 4;
  geodesic_draping::GeoDrapeSolver refineSolver(mesh, fixtureHeatOptions(), {}, refine);
  requireFiniteFaceShear(refineSolver.solve(seedXY, angleDegrees, options), "delaunay refine smoke");
}

void testIntrinsicRetrieval(const std::filesystem::path& root) {
  const std::filesystem::path fixtureDir = root / "small_curved";
  const auto mesh = geodesic_draping::fixture_io::loadMesh(fixtureDir);
  const auto seedXY = geodesic_draping::fixture_io::loadSeedXY(fixtureDir);
  const double angleDegrees = geodesic_draping::fixture_io::loadAngleDegrees(fixtureDir);

  geodesic_draping::DrapeSolveOptions options;
  options.mode = geodesic_draping::DrapeSolveMode::Hybrid;
  options.retrieval = geodesic_draping::ResultDomain::Intrinsic;
  options.sampleVertexShear = true;

  geodesic_draping::GeoDrapeSolver solver(mesh, fixtureHeatOptions());
  const auto result = solver.solve(seedXY, angleDegrees, options);
  assert(result.domain == geodesic_draping::ResultDomain::Intrinsic);
  assert(result.mesh.domain == geodesic_draping::ResultDomain::Intrinsic);
  assert(!result.mesh.vertices3D);
  assert(result.mesh.edgeLengths);
  assert(result.mesh.gluingMap);
  assert(result.mesh.faces.size() == result.mesh.edgeLengths->size());
  assert(result.mesh.faces.size() == result.mesh.gluingMap->size());
  assert(result.origin.intrinsicPoint);
  assert(result.origin.intrinsicFamilyDirections[0]);
  assert(result.origin.intrinsicFamilyDirections[1]);
  assert(!result.traces[0].positive.intrinsicPoints.empty());
  assert(!result.traces[0].negative.intrinsicPoints.empty());
  assert(result.distances);
  assert(result.faceShear);
  assert(result.vertexShear);
}

void testSubdivisionRetrieval(const std::filesystem::path& root) {
  const std::filesystem::path fixtureDir = root / "small_curved";
  const auto mesh = geodesic_draping::fixture_io::loadMesh(fixtureDir);
  const auto seedXY = geodesic_draping::fixture_io::loadSeedXY(fixtureDir);
  const double angleDegrees = geodesic_draping::fixture_io::loadAngleDegrees(fixtureDir);

  geodesic_draping::RefinementOptions refinement;
  refinement.mode = geodesic_draping::RefinementMode::DelaunayFlip;

  geodesic_draping::DrapeSolveOptions options;
  options.mode = geodesic_draping::DrapeSolveMode::Hybrid;
  options.retrieval = geodesic_draping::ResultDomain::Subdivision;
  options.sampleVertexShear = true;

  geodesic_draping::GeoDrapeSolver solver(mesh, fixtureHeatOptions(), {}, refinement);
  const auto result = solver.solve(seedXY, angleDegrees, options);
  assert(result.domain == geodesic_draping::ResultDomain::Subdivision);
  assert(result.mesh.domain == geodesic_draping::ResultDomain::Subdivision);
  assert(result.mesh.vertices3D);
  assert(!result.mesh.vertices3D->empty());
  assert(!result.mesh.faces.empty());
  assert(result.faceShear);
  assert(result.faceShear->size() == result.mesh.faces.size());
  assert(result.vertexShear);
  assert(result.vertexShear->size() == result.mesh.vertices3D->size());
  assert(result.distances);
  assert((*result.distances)[0].size() == result.mesh.vertices3D->size());
  assert(result.directions[0].size() == result.mesh.faces.size());
  assert(!result.traces[0].positive.extrinsicPoints.empty());
  for (double value : *result.faceShear) {
    assert(std::isfinite(value));
  }
}

} // namespace

int main() {
  const std::filesystem::path fixtureRoot = GEODESIC_DRAPING_TEST_DATA_DIR;

  const auto tiny = solveFixture(fixtureRoot, "tiny_planar");
  requireFiniteComplete(tiny, geodesic_draping::fixture_io::loadMesh(fixtureRoot / "tiny_planar").vertices.size());
  requireNearArray((*tiny.distances)[0],
                   geodesic_draping::fixture_io::loadGoldenScalarArray(fixtureRoot / "tiny_planar", "dist_0"),
                   1e-8,
                   "tiny complete dist_0");
  requireNearArray((*tiny.distances)[1],
                   geodesic_draping::fixture_io::loadGoldenScalarArray(fixtureRoot / "tiny_planar", "dist_1"),
                   1e-8,
                   "tiny complete dist_1");
  requireNearArray(*tiny.vertexShear,
                   geodesic_draping::fixture_io::loadGoldenShearArray(fixtureRoot / "tiny_planar", "complete"),
                   1e-8,
                   "tiny complete shear");

  const auto small = solveFixture(fixtureRoot, "small_curved");
  requireFiniteComplete(small, geodesic_draping::fixture_io::loadMesh(fixtureRoot / "small_curved").vertices.size());

  const auto smoothGood = solveFixture(fixtureRoot, "smooth_quality_good");
  requireFiniteComplete(smoothGood,
                geodesic_draping::fixture_io::loadMesh(fixtureRoot / "smooth_quality_good").vertices.size());

  const auto smoothPoor = solveFixture(fixtureRoot, "smooth_quality_poor");
  requireFiniteComplete(smoothPoor,
                geodesic_draping::fixture_io::loadMesh(fixtureRoot / "smooth_quality_poor").vertices.size());

  const auto demo = solveFixture(fixtureRoot, "demo_part");
  requireFiniteComplete(demo, geodesic_draping::fixture_io::loadMesh(fixtureRoot / "demo_part").vertices.size());
  testHighLevelGeneratorGoldenEndpoints(fixtureRoot, "demo_part");
  testPersistentSolver(fixtureRoot, "demo_part");
  testOneShotComplete(fixtureRoot, "tiny_planar");
  try {
    testRefinementSmoke(fixtureRoot);
    testIntrinsicRetrieval(fixtureRoot);
    testSubdivisionRetrieval(fixtureRoot);
  } catch (const std::exception& e) {
    std::cerr << "architecture smoke failed: " << e.what() << "\n";
    return 1;
  }

  return 0;
}

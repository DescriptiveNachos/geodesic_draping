#include "fixture_io.h"
#include "geodesic_draping/geodrape.h"
#include "geodesic_draping/quality.h"

#include <cassert>
#include <cmath>
#include <fstream>
#include <filesystem>
#include <stdexcept>

namespace {

geodesic_draping::SignedHeatSolveOptions fixtureHeatOptions() {
  geodesic_draping::SignedHeatSolveOptions options;
  options.preserveSourceNormals = false;
  options.levelSetConstraint = geometrycentral::LevelSetConstraint::None;
  options.softLevelSetWeight = -1.0;
  return options;
}

geodesic_draping::DrapeResult solveFixture(const std::filesystem::path& root,
                                           const std::string& name,
                                           geodesic_draping::DrapeSolveMode mode) {
  const std::filesystem::path fixtureDir = root / name;
  geodesic_draping::DrapeSolveOptions options;
  options.mode = mode;
  return geodesic_draping::solveDrape(
      geodesic_draping::fixture_io::loadMesh(fixtureDir),
      geodesic_draping::fixture_io::loadSeedXY(fixtureDir),
      geodesic_draping::fixture_io::loadAngleDegrees(fixtureDir),
      fixtureHeatOptions(),
      options);
}

void testFixtureReports(const std::filesystem::path& root) {
  const std::filesystem::path tinyDir = root / "tiny_planar";
  const auto tinyMesh = geodesic_draping::fixture_io::loadMesh(tinyDir);
  const auto tiny = solveFixture(root, "tiny_planar", geodesic_draping::DrapeSolveMode::Complete);
  const auto report = geodesic_draping::analyzeSolveQuality(tinyMesh, tiny);

  assert(report.mesh.edgeLengthP99ToMedian >= 1.0);
  assert(report.mesh.triangleAspectRatioMax > 0.0);
  assert(report.generators.lengths.size() == 4);
  assert(report.primaryShear.sampleCount == tinyMesh.vertices.size());
  assert(report.primaryShear.nonFiniteCount == 0);
}

void testAspectRatioWarning() {
  geodesic_draping::SurfaceMeshData mesh;
  mesh.vertices = {
      {0.0, 0.0, 0.0},
      {100.0, 0.0, 0.0},
      {0.001, 0.001, 0.0},
      {0.0, 1.0, 0.0},
  };
  mesh.faces = {{{0, 1, 2}}, {{0, 2, 3}}};

  geodesic_draping::DrapeResult result;
  result.mode = geodesic_draping::DrapeSolveMode::Fast;
  result.origin.extrinsicPoint = geodesic_draping::Vec3{0.1, 0.1, 0.0};
  for (auto& family : result.traces) {
    family.positive.hitBoundary = true;
    family.positive.length = 1.0;
    family.negative.hitBoundary = true;
    family.negative.length = 1.0;
  }
  result.faceShearAnglesDegrees = std::vector<double>{0.0, 0.0};

  const auto report = geodesic_draping::analyzeSolveQuality(mesh, result);
  assert(report.level != geodesic_draping::SolveQualityLevel::Good);
  assert(report.mesh.triangleAspectRatioMax > 1000.0);
}

void testDemoPartLocalShearWarning(const std::filesystem::path& root) {
  const std::filesystem::path fixtureDir = root / "demo_part";
  const auto mesh = geodesic_draping::fixture_io::loadMesh(fixtureDir);
  const auto result = solveFixture(root, "demo_part", geodesic_draping::DrapeSolveMode::Fast);
  auto thresholds = geodesic_draping::SolveQualityThresholds{};
  thresholds.localShearJumpMaxWarning = 40.0;
  thresholds.localShearJumpMaxPoor = 1000.0;
  thresholds.localShearJumpP99Warning = 1000.0;
  thresholds.shearMaxMinusP99Warning = 1000.0;
  const auto report = geodesic_draping::analyzeSolveQuality(mesh, result, thresholds);
  assert(report.primaryShear.localJumpMax > thresholds.localShearJumpMaxWarning);
  assert(report.level == geodesic_draping::SolveQualityLevel::Warning);
}

void testThresholdConfigLoader() {
  const std::filesystem::path path = std::filesystem::current_path() / "quality_thresholds_test.cfg";
  {
    std::ofstream output(path);
    output << "# partial quality threshold override\n";
    output << "triangleAspectRatioWarning = 12.5\n";
    output << "localShearJumpMaxPoor = 90\n";
  }

  const auto thresholds = geodesic_draping::loadSolveQualityThresholds(path);
  assert(thresholds.triangleAspectRatioWarning == 12.5);
  assert(thresholds.localShearJumpMaxPoor == 90.0);
  assert(thresholds.edgeLengthP99ToMedianWarning == geodesic_draping::SolveQualityThresholds{}.edgeLengthP99ToMedianWarning);
  std::filesystem::remove(path);

  const std::filesystem::path badPath = std::filesystem::current_path() / "quality_thresholds_bad_test.cfg";
  {
    std::ofstream output(badPath);
    output << "notAThreshold = 1\n";
  }

  bool threw = false;
  try {
    (void)geodesic_draping::loadSolveQualityThresholds(badPath);
  } catch (const std::runtime_error&) {
    threw = true;
  }
  std::filesystem::remove(badPath);
  assert(threw);
}

} // namespace

int main() {
  const std::filesystem::path fixtureRoot = GEODESIC_DRAPING_TEST_DATA_DIR;
  testFixtureReports(fixtureRoot);
  testAspectRatioWarning();
  testDemoPartLocalShearWarning(fixtureRoot);
  testThresholdConfigLoader();
  return 0;
}

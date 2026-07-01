#include "fixture_io.h"
#include "geodesic_draping/generator_tracing.h"
#include "geodesic_draping/seed_projection.h"

#include <cassert>
#include <cmath>
#include <filesystem>
#include <iostream>
#include <string>

namespace {

double maxAbsDiff(const geodesic_draping::Vec3& a, const geodesic_draping::Vec3& b) {
  return (a - b).cwiseAbs().maxCoeff();
}

void requireNear(const geodesic_draping::Vec3& a,
                 const geodesic_draping::Vec3& b,
                 double tolerance,
                 const std::string& label) {
  if (maxAbsDiff(a, b) > tolerance) {
    std::cerr << label << " mismatch\n"
              << "actual:   " << a.transpose() << "\n"
              << "expected: " << b.transpose() << "\n"
              << "diff:     " << maxAbsDiff(a, b) << "\n";
    assert(false);
  }
}

void testFixture(const std::filesystem::path& root, const std::string& name, bool compareGoldenEnds) {
  const std::filesystem::path fixtureDir = root / name;
  const geodesic_draping::SurfaceMeshData meshData = geodesic_draping::fixture_io::loadMesh(fixtureDir);
  auto surface = geodesic_draping::makeGeometryCentralSurface(meshData);

  const geodesic_draping::Vec2 seedXY = geodesic_draping::fixture_io::loadSeedXY(fixtureDir);
  const double angleDegrees = geodesic_draping::fixture_io::loadAngleDegrees(fixtureDir);
  const auto projection = geodesic_draping::projectPointXYToMesh(meshData, seedXY);
  assert(projection.has_value());
  const auto directions = geodesic_draping::generateOrthogonalDirections(angleDegrees);

  const auto traces = geodesic_draping::traceGenerators(surface, projection->surfacePoint, directions);
  const auto sourceCurves = geodesic_draping::pairOppositeGeneratorTraces(traces);
  const auto goldenEnds = geodesic_draping::fixture_io::loadGoldenGeneratorLastPoints(fixtureDir);
  const auto goldenPairedCounts = geodesic_draping::fixture_io::loadGoldenPairedGeneratorPointCounts(fixtureDir);
  assert(goldenEnds.size() == traces.size());

  for (size_t i = 0; i < traces.size(); ++i) {
    assert(!traces[i].points.empty());
    assert(traces[i].points.size() == traces[i].surfaceReferences.size());
    assert(traces[i].length > 0.0);
    requireNear(traces[i].points.front(), projection->cartesian, 1e-9,
                name + " trace " + std::to_string(i) + " start");

    if (compareGoldenEnds) {
      requireNear(traces[i].points.back(), goldenEnds[i], 1e-8,
                  name + " trace " + std::to_string(i) + " end");
    }
  }

  assert(sourceCurves.curves[0].size() == traces[1].surfaceReferences.size() + traces[0].surfaceReferences.size());
  assert(sourceCurves.curves[1].size() == traces[3].surfaceReferences.size() + traces[2].surfaceReferences.size());
  if (compareGoldenEnds) {
    assert(sourceCurves.curves[0].size() == goldenPairedCounts[0]);
    assert(sourceCurves.curves[1].size() == goldenPairedCounts[1]);
  }
}

} // namespace

int main() {
  const std::filesystem::path fixtureRoot = GEODESIC_DRAPING_TEST_DATA_DIR;
  testFixture(fixtureRoot, "tiny_planar", true);
  testFixture(fixtureRoot, "small_curved", false);
  testFixture(fixtureRoot, "demo_part", false);
  return 0;
}

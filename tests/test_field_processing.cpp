#include "fixture_io.h"
#include "geodesic_draping/field_processing.h"

#include <cassert>
#include <cmath>
#include <filesystem>
#include <iostream>
#include <string>
#include <vector>

namespace {

void requireNearArray(const std::vector<double>& actual,
                      const std::vector<double>& golden,
                      double tolerance,
                      const std::string& label) {
  assert(actual.size() == golden.size());
  double maxDiff = 0.0;
  for (size_t i = 0; i < actual.size(); ++i) {
    if (std::isfinite(actual[i]) || std::isfinite(golden[i])) {
      maxDiff = std::max(maxDiff, std::abs(actual[i] - golden[i]));
    }
  }
  if (maxDiff > tolerance) {
    std::cerr << label << " max diff " << maxDiff << " exceeds tolerance " << tolerance << "\n";
    assert(false);
  }
}

void testFixture(const std::filesystem::path& root, const std::string& name) {
  const std::filesystem::path fixtureDir = root / name;

  const auto completeShear = geodesic_draping::computeShearAnglesDegrees(
      geodesic_draping::fixture_io::loadGoldenVectorArray(fixtureDir, "complete", "grad_0"),
      geodesic_draping::fixture_io::loadGoldenVectorArray(fixtureDir, "complete", "grad_1"));
  requireNearArray(completeShear,
                   geodesic_draping::fixture_io::loadGoldenShearArray(fixtureDir, "complete"),
                   1e-10,
                   name + " complete shear");

  const auto fastShear = geodesic_draping::computeShearAnglesDegrees(
      geodesic_draping::fixture_io::loadGoldenVectorArray(fixtureDir, "fast", "grad_0"),
      geodesic_draping::fixture_io::loadGoldenVectorArray(fixtureDir, "fast", "grad_1"));
  requireNearArray(fastShear,
                   geodesic_draping::fixture_io::loadGoldenShearArray(fixtureDir, "fast"),
                   1e-10,
                   name + " fast shear");
}

} // namespace

int main() {
  const std::filesystem::path fixtureRoot = GEODESIC_DRAPING_TEST_DATA_DIR;
  testFixture(fixtureRoot, "tiny_planar");
  testFixture(fixtureRoot, "small_curved");
  testFixture(fixtureRoot, "demo_part");
  return 0;
}

#include "fixture_io.h"
#include "geodesic_draping/field_processing.h"

#include <cassert>
#include <cmath>
#include <filesystem>
#include <iostream>
#include <stdexcept>
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

void testGradientFixture(const std::filesystem::path& root,
                         const std::string& name,
                         double tolerance) {
  const std::filesystem::path fixtureDir = root / name;
  const auto mesh = geodesic_draping::fixture_io::loadMesh(fixtureDir);
  const auto gradients = geodesic_draping::computeVertexScalarGradients(
      mesh,
      geodesic_draping::fixture_io::loadGoldenScalarArray(fixtureDir, "dist_0"));
  requireNearVectorArray(gradients,
                         geodesic_draping::fixture_io::loadGoldenVectorArray(fixtureDir, "complete", "grad_0"),
                         tolerance,
                         name + " grad_0");
}

void testFaceScalarAveraging() {
  geodesic_draping::SurfaceMeshData mesh;
  mesh.vertices = {
      {0.0, 0.0, 0.0},
      {1.0, 0.0, 0.0},
      {1.0, 1.0, 0.0},
      {0.0, 1.0, 0.0},
  };
  mesh.faces = {{{0, 1, 2}}, {{0, 2, 3}}};

  const std::vector<double> averaged = geodesic_draping::averageFaceScalarsToVertices(
      mesh,
      {2.0, 4.0},
      geodesic_draping::FaceScalarAveraging::FaceArea);
  requireNearArray(averaged, {3.0, 2.0, 3.0, 4.0}, 1e-12, "area face scalar average");

  bool threw = false;
  try {
    (void)geodesic_draping::averageFaceScalarsToVertices(mesh, {1.0});
  } catch (const std::runtime_error&) {
    threw = true;
  }
  assert(threw);
}

void testVertexShearAngles() {
  using geodesic_draping::Vec3;
  const auto orthogonal = geodesic_draping::computeShearAnglesDegrees(
      {Vec3{1.0, 0.0, 0.0}},
      {Vec3{0.0, 1.0, 0.0}});
  requireNearArray(orthogonal, {0.0}, 1e-12, "orthogonal vertex shear");

  const auto parallel = geodesic_draping::computeShearAnglesDegrees(
      {Vec3{1.0, 0.0, 0.0}},
      {Vec3{2.0, 0.0, 0.0}});
  requireNearArray(parallel, {90.0}, 1e-12, "parallel vertex shear");

  const auto opposite = geodesic_draping::computeShearAnglesDegrees(
      {Vec3{1.0, 0.0, 0.0}},
      {Vec3{-2.0, 0.0, 0.0}});
  requireNearArray(opposite, {90.0}, 1e-12, "opposite vertex shear");
}

} // namespace

int main() {
  const std::filesystem::path fixtureRoot = GEODESIC_DRAPING_TEST_DATA_DIR;
  testFaceScalarAveraging();
  testVertexShearAngles();
  testFixture(fixtureRoot, "tiny_planar");
  testFixture(fixtureRoot, "small_curved");
  testFixture(fixtureRoot, "demo_part");
  testGradientFixture(fixtureRoot, "tiny_planar", 1e-12);
  testGradientFixture(fixtureRoot, "small_curved", 2e-3);
  testGradientFixture(fixtureRoot, "demo_part", 1e-1);
  return 0;
}

#include "fixture_io.h"
#include "geodesic_draping/geodrape.h"
#include "geodesic_draping/plotting.h"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

using geodesic_draping::CompleteDrapeResult;
using geodesic_draping::FastDrapeResult;
using geodesic_draping::ProjectionPlotOptions;
using geodesic_draping::Vec3;

void printUsage(const char* argv0) {
  std::cout << "Usage:\n"
            << "  " << argv0 << " [fixture-name] [--no-show] [--direction-length value]\n\n"
            << "Fixtures are loaded from " << GEODESIC_DRAPING_TEST_DATA_DIR << "\n"
            << "Default fixture: demo_part\n";
}

double maxAbsDiff(const Vec3& a, const Vec3& b) {
  return (a - b).cwiseAbs().maxCoeff();
}

double maxAbsDiff(const std::vector<double>& actual, const std::vector<double>& golden) {
  const size_t n = std::min(actual.size(), golden.size());
  double diff = 0.0;
  for (size_t i = 0; i < n; ++i) {
    diff = std::max(diff, std::abs(actual[i] - golden[i]));
  }
  return diff;
}

struct ScalarStats {
  double min = 0.0;
  double max = 0.0;
  double mean = 0.0;
};

ScalarStats stats(const std::vector<double>& values) {
  if (values.empty()) {
    return {};
  }
  ScalarStats s;
  s.min = std::numeric_limits<double>::infinity();
  s.max = -std::numeric_limits<double>::infinity();
  double sum = 0.0;
  size_t count = 0;
  for (double value : values) {
    if (!std::isfinite(value)) {
      continue;
    }
    s.min = std::min(s.min, value);
    s.max = std::max(s.max, value);
    sum += value;
    ++count;
  }
  s.mean = count > 0 ? sum / static_cast<double>(count) : std::numeric_limits<double>::quiet_NaN();
  return s;
}

void printStats(const std::string& label, const std::vector<double>& values) {
  const ScalarStats s = stats(values);
  std::cout << std::left << std::setw(18) << label
            << " min " << s.min
            << "  max " << s.max
            << "  mean " << s.mean << "\n";
}

void printVectorComparison(const std::string& label, const Vec3& actual, const Vec3& golden) {
  std::cout << std::left << std::setw(18) << label << " actual [" << actual.transpose() << "]"
            << "  golden [" << golden.transpose() << "]"
            << "  max_abs_diff " << maxAbsDiff(actual, golden) << "\n";
}

void printScalarComparison(const std::string& label,
                           const std::vector<double>& actual,
                           const std::vector<double>& golden) {
  std::cout << std::left << std::setw(18) << label
            << " actual_n " << actual.size()
            << "  golden_n " << golden.size()
            << "  max_abs_diff " << maxAbsDiff(actual, golden) << "\n";
}

void printVectorArrayComparison(const std::string& label,
                                const std::vector<Vec3>& actual,
                                const std::vector<Vec3>& golden) {
  const size_t n = std::min(actual.size(), golden.size());
  double diff = 0.0;
  for (size_t i = 0; i < n; ++i) {
    diff = std::max(diff, maxAbsDiff(actual[i], golden[i]));
  }
  std::cout << std::left << std::setw(18) << label
            << " actual_n " << actual.size()
            << "  golden_n " << golden.size()
            << "  max_abs_diff " << diff << "\n";
}

geodesic_draping::SignedHeatSolveOptions fixtureHeatOptions() {
  geodesic_draping::SignedHeatSolveOptions options;
  options.preserveSourceNormals = false;
  options.levelSetConstraint = geometrycentral::LevelSetConstraint::None;
  options.softLevelSetWeight = -1.0;
  return options;
}

} // namespace

int main(int argc, char** argv) {
  try {
    std::string fixtureName = "demo_part";
    bool show = true;
    double directionLength = 25.0;

    for (int i = 1; i < argc; ++i) {
      const std::string arg = argv[i];
      if (arg == "--help" || arg == "-h") {
        printUsage(argv[0]);
        return 0;
      }
      if (arg == "--no-show") {
        show = false;
        continue;
      }
      if (arg == "--direction-length") {
        if (i + 1 >= argc) {
          throw std::runtime_error("--direction-length requires a value");
        }
        directionLength = std::stod(argv[++i]);
        continue;
      }
      fixtureName = arg;
    }

    const std::filesystem::path fixtureRoot = GEODESIC_DRAPING_TEST_DATA_DIR;
    const std::filesystem::path fixtureDir = fixtureRoot / fixtureName;
    const auto mesh = geodesic_draping::fixture_io::loadMesh(fixtureDir);
    const auto seedXY = geodesic_draping::fixture_io::loadSeedXY(fixtureDir);
    const double angleDegrees = geodesic_draping::fixture_io::loadAngleDegrees(fixtureDir);
    const auto heatOptions = fixtureHeatOptions();
    const CompleteDrapeResult result =
        geodesic_draping::solveCompleteDrape(mesh, seedXY, angleDegrees, heatOptions);
    const FastDrapeResult fastResult =
        geodesic_draping::solveFastDrape(mesh, seedXY, angleDegrees, heatOptions);

    const Vec3 goldenOrigin = geodesic_draping::fixture_io::loadGoldenOrigin(fixtureDir);
    const size_t goldenFaceIndex = geodesic_draping::fixture_io::loadGoldenSeedFaceIndex(fixtureDir);
    const Vec3 goldenBarycentric = geodesic_draping::fixture_io::loadGoldenSeedBarycentric(fixtureDir);
    const std::vector<Vec3> goldenDirections = geodesic_draping::fixture_io::loadGoldenDirections(fixtureDir);
    const std::vector<Vec3> goldenGeneratorEnds =
        geodesic_draping::fixture_io::loadGoldenGeneratorLastPoints(fixtureDir);
    const std::vector<size_t> goldenGeneratorCounts =
        geodesic_draping::fixture_io::loadGoldenGeneratorPointCounts(fixtureDir);
    const std::vector<size_t> goldenPairedCounts =
        geodesic_draping::fixture_io::loadGoldenPairedGeneratorPointCounts(fixtureDir);

    std::cout << std::setprecision(17);
    std::cout << "Fixture: " << fixtureName << "\n"
              << "Vertices: " << mesh.vertices.size() << "  Faces: " << mesh.faces.size() << "\n"
              << "Seed XY: [" << seedXY.transpose() << "]  Angle degrees: " << angleDegrees << "\n\n";

    printVectorComparison("origin", result.seed.cartesian, goldenOrigin);
    std::cout << std::left << std::setw(18) << "face_index"
              << " actual " << result.seed.surfacePoint.faceIndex
              << "  golden " << goldenFaceIndex;
    if (result.seed.surfacePoint.faceIndex != goldenFaceIndex) {
      std::cout << "  note: valid tie-order differences can occur on edges/vertices";
    }
    std::cout << "\n";
    printVectorComparison("barycentric", result.seed.surfacePoint.barycentric, goldenBarycentric);
    std::cout << "\n";

    for (size_t i = 0; i < result.directions.size(); ++i) {
      printVectorComparison("direction " + std::to_string(i), result.directions[i], goldenDirections[i]);
    }
    std::cout << "\n";

    for (size_t i = 0; i < result.generators.size(); ++i) {
      std::cout << "generator " << i
                << " actual_points " << result.generators[i].points.size()
                << "  golden_points " << goldenGeneratorCounts[i]
                << "  hit_boundary " << (result.generators[i].hitBoundary ? "true" : "false")
                << "  length " << result.generators[i].length << "\n";
      printVectorComparison("  end", result.generators[i].points.back(), goldenGeneratorEnds[i]);
    }
    std::cout << "\n";

    for (size_t i = 0; i < result.sourceCurves.curves.size(); ++i) {
      std::cout << "source curve " << i
                << " actual_refs " << result.sourceCurves.curves[i].size()
                << "  golden_refs " << goldenPairedCounts[i] << "\n";
    }
    std::cout << "\n";

    printStats("dist_0", result.distances[0]);
    printStats("dist_1", result.distances[1]);
    printStats("shear_degrees", result.shearAnglesDegrees);
    printScalarComparison("dist_0",
                          result.distances[0],
                          geodesic_draping::fixture_io::loadGoldenScalarArray(fixtureDir, "dist_0"));
    printScalarComparison("dist_1",
                          result.distances[1],
                          geodesic_draping::fixture_io::loadGoldenScalarArray(fixtureDir, "dist_1"));
    printScalarComparison("shear",
                          result.shearAnglesDegrees,
                          geodesic_draping::fixture_io::loadGoldenShearArray(fixtureDir, "complete"));
    std::cout << "\n";
    printStats("fast_shear", fastResult.shearAnglesDegrees);
    printVectorArrayComparison("fast grad_0",
                               fastResult.gradients[0],
                               geodesic_draping::fixture_io::loadGoldenVectorArray(fixtureDir, "fast", "grad_0"));
    printVectorArrayComparison("fast grad_1",
                               fastResult.gradients[1],
                               geodesic_draping::fixture_io::loadGoldenVectorArray(fixtureDir, "fast", "grad_1"));
    printScalarComparison("fast shear",
                          fastResult.shearAnglesDegrees,
                          geodesic_draping::fixture_io::loadGoldenShearArray(fixtureDir, "fast"));

    ProjectionPlotOptions options;
    options.name = fixtureName + " drape comparison";
    options.directionLength = directionLength;
    options.clearExisting = true;
    options.show = show;
    geodesic_draping::plotDrapeComparisonResult(mesh, result, fastResult, options);

    if (!show) {
      std::cout << "\nPlot data registered with Polyscope, but --no-show was set.\n";
    }
    return 0;
  } catch (const std::exception& e) {
    std::cerr << "error: " << e.what() << "\n\n";
    printUsage(argv[0]);
    return 1;
  }
}

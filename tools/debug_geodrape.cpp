#include "fixture_io.h"
#include "geodesic_draping/diagnostics.h"
#include "geodesic_draping/field_processing.h"
#include "geodesic_draping/geodrape.h"
#include "geodesic_draping/plotting.h"
#include "geodesic_draping/quality.h"

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

using geodesic_draping::MagnitudeStats;
using geodesic_draping::ProjectionPlotOptions;
using geodesic_draping::SolveQualityReport;
using geodesic_draping::Vec3;
using geodesic_draping::VectorMagnitudeDiagnostics;

void printUsage(const char* argv0) {
  std::cout << "Usage:\n"
            << "  " << argv0
            << " [fixture-name] [--no-show] [--direction-length value] [--quality-thresholds path]\n\n"
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

void printMagnitudeDiagnostics(const std::string& label,
                               const VectorMagnitudeDiagnostics& diagnostics) {
  const MagnitudeStats& s = diagnostics.stats;
  std::cout << std::left << std::setw(18) << label
            << " |v| min " << s.min
            << "  max " << s.max
            << "  mean " << s.mean
            << "  max_abs_unit_dev " << s.maxAbsDeviationFromUnit
            << "  mean_abs_unit_dev " << s.meanAbsDeviationFromUnit
            << "  near_zero " << s.nearZeroCount
            << "  nonfinite " << s.nonFiniteCount << "\n";
}

void printShearSummary(const std::string& label, const std::vector<double>& shear) {
  const ScalarStats s = stats(shear);
  std::cout << std::left << std::setw(18) << label
            << " shear min " << s.min
            << "  max " << s.max
            << "  mean " << s.mean << "\n";
}

void printQualityReport(const std::string& label, const SolveQualityReport& report) {
  std::cout << "\n" << label << " quality " << geodesic_draping::toString(report.level) << "\n";
  std::cout << "  mesh edge_p99/median " << report.mesh.edgeLengthP99ToMedian
            << "  aspect_p99 " << report.mesh.triangleAspectRatioP99
            << "  aspect_max " << report.mesh.triangleAspectRatioMax
            << "  min_area/bbox " << report.mesh.minFaceAreaRelativeToBbox << "\n";
  std::cout << "  generators nearest_wall " << report.generators.nearestWallDistance
            << "  farthest_wall " << report.generators.farthestWallDistance
            << "  plausible_length [" << report.generators.plausibleMinLength
            << ", " << report.generators.plausibleMaxLength << "]"
            << "  short " << report.generators.shortGeneratorCount
            << "  long " << report.generators.longGeneratorCount
            << "  missed_boundary " << report.generators.missedBoundaryCount << "\n";
  std::cout << "  primary shear p95 " << report.primaryShear.p95
            << "  p99 " << report.primaryShear.p99
            << "  max-p99 " << report.primaryShear.maxMinusP99
            << "  local_jump_p99 " << report.primaryShear.localJumpP99
            << "  local_jump_max " << report.primaryShear.localJumpMax
            << "  nonfinite " << report.primaryShear.nonFiniteCount << "\n";
  for (const std::string& warning : report.warnings) {
    std::cout << "  warning: " << warning << "\n";
  }
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
    geodesic_draping::SolveQualityThresholds qualityThresholds;

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
      if (arg == "--quality-thresholds") {
        if (i + 1 >= argc) {
          throw std::runtime_error("--quality-thresholds requires a path");
        }
        qualityThresholds = geodesic_draping::loadSolveQualityThresholds(argv[++i], qualityThresholds);
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
    geodesic_draping::DrapeSolveOptions completeOptions;
    completeOptions.mode = geodesic_draping::DrapeSolveMode::Complete;
    const geodesic_draping::DrapeResult result =
        geodesic_draping::solveDrape(mesh, seedXY, angleDegrees, heatOptions, completeOptions);
    geodesic_draping::DrapeSolveOptions completeFieldOptions = completeOptions;
    completeFieldOptions.retrieval = geodesic_draping::ResultDomain::Intrinsic;
    const geodesic_draping::DrapeResult completeFieldResult =
        geodesic_draping::solveDrape(mesh, seedXY, angleDegrees, heatOptions, completeFieldOptions);
    geodesic_draping::DrapeSolveOptions fastOptions;
    fastOptions.mode = geodesic_draping::DrapeSolveMode::Fast;
    fastOptions.retrieval = geodesic_draping::ResultDomain::Intrinsic;
    fastOptions.sampleVertexShear = true;
    const geodesic_draping::DrapeResult fastResult =
        geodesic_draping::solveDrape(mesh, seedXY, angleDegrees, heatOptions, fastOptions);
    if (!result.distances || !result.vertexShear ||
        !completeFieldResult.distances || !completeFieldResult.vertexShear ||
        !fastResult.faceShear || !fastResult.vertexShear) {
      throw std::runtime_error("debug solve did not return the expected fields");
    }
    if (!result.origin.intrinsicPoint ||
        result.origin.intrinsicPoint->params.size() != 3 ||
        !result.origin.extrinsicPoint ||
        !result.origin.extrinsicFamilyDirections[0] ||
        !result.origin.extrinsicFamilyDirections[1]) {
      throw std::runtime_error("debug solve did not return the expected origin and direction references");
    }

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

    printVectorComparison("origin", *result.origin.extrinsicPoint, goldenOrigin);
    std::cout << std::left << std::setw(18) << "face_index"
              << " actual " << result.origin.intrinsicPoint->elementIndex
              << "  golden " << goldenFaceIndex;
    if (result.origin.intrinsicPoint->elementIndex != goldenFaceIndex) {
      std::cout << "  note: valid tie-order differences can occur on edges/vertices";
    }
    std::cout << "\n";
    printVectorComparison("barycentric",
                          Vec3(result.origin.intrinsicPoint->params[0],
                               result.origin.intrinsicPoint->params[1],
                               result.origin.intrinsicPoint->params[2]),
                          goldenBarycentric);
    std::cout << "\n";

    const std::array<Vec3, 4> retrievedDirections = {
        *result.origin.extrinsicFamilyDirections[0],
        -*result.origin.extrinsicFamilyDirections[0],
        *result.origin.extrinsicFamilyDirections[1],
        -*result.origin.extrinsicFamilyDirections[1],
    };
    for (size_t i = 0; i < retrievedDirections.size(); ++i) {
      printVectorComparison("direction " + std::to_string(i), retrievedDirections[i], goldenDirections[i]);
    }
    std::cout << "\n";

    const std::array<const geodesic_draping::DrapeTrace*, 4> retrievedTraces = {
        &result.traces[0].positive,
        &result.traces[0].negative,
        &result.traces[1].positive,
        &result.traces[1].negative,
    };
    for (size_t i = 0; i < retrievedTraces.size(); ++i) {
      std::cout << "generator " << i
                << " actual_points " << retrievedTraces[i]->extrinsicPoints.size()
                << "  golden_points " << goldenGeneratorCounts[i]
                << "  hit_boundary " << (retrievedTraces[i]->hitBoundary ? "true" : "false")
                << "  length " << retrievedTraces[i]->length << "\n";
      printVectorComparison("  end", retrievedTraces[i]->extrinsicPoints.back(), goldenGeneratorEnds[i]);
    }
    std::cout << "\n";

    for (size_t i = 0; i < result.traces.size(); ++i) {
      const size_t actualRefs =
          result.traces[i].negative.extrinsicPoints.size() + result.traces[i].positive.extrinsicPoints.size();
      std::cout << "source curve " << i
                << " actual_refs " << actualRefs
                << "  golden_refs " << goldenPairedCounts[i] << "\n";
    }
    std::cout << "\n";

    printStats("dist_0", (*result.distances)[0]);
    printStats("dist_1", (*result.distances)[1]);
    printStats("shear_degrees", *result.vertexShear);
    printScalarComparison("dist_0",
                          (*result.distances)[0],
                          geodesic_draping::fixture_io::loadGoldenScalarArray(fixtureDir, "dist_0"));
    printScalarComparison("dist_1",
                          (*result.distances)[1],
                          geodesic_draping::fixture_io::loadGoldenScalarArray(fixtureDir, "dist_1"));
    printScalarComparison("shear",
                          *result.vertexShear,
                          geodesic_draping::fixture_io::loadGoldenShearArray(fixtureDir, "complete"));
    std::cout << "\n";
    printStats("fast_face_shear", *fastResult.faceShear);
    printStats("fast_vertex_shear", *fastResult.vertexShear);
    std::cout << "\n";
    auto surface = geodesic_draping::makeGeometryCentralSurface(mesh);
    const auto completeDirections0 = geodesic_draping::directionsToExtrinsicVectors(surface, completeFieldResult.directions[0]);
    const auto completeDirections1 = geodesic_draping::directionsToExtrinsicVectors(surface, completeFieldResult.directions[1]);
    printMagnitudeDiagnostics("complete face dir_0", geodesic_draping::analyzeVectorMagnitudes(completeDirections0));
    printMagnitudeDiagnostics("complete face dir_1", geodesic_draping::analyzeVectorMagnitudes(completeDirections1));
    const auto directions0 = geodesic_draping::directionsToExtrinsicVectors(surface, fastResult.directions[0]);
    const auto directions1 = geodesic_draping::directionsToExtrinsicVectors(surface, fastResult.directions[1]);
    printMagnitudeDiagnostics("fast face dir_0", geodesic_draping::analyzeVectorMagnitudes(directions0));
    printMagnitudeDiagnostics("fast face dir_1", geodesic_draping::analyzeVectorMagnitudes(directions1));
    printShearSummary("fast face", *fastResult.faceShear);
    printShearSummary("fast vertex", *fastResult.vertexShear);
    printQualityReport("complete", geodesic_draping::analyzeSolveQuality(mesh, result, qualityThresholds));

    ProjectionPlotOptions options;
    options.name = fixtureName + " drape comparison";
    options.directionLength = directionLength;
    options.clearExisting = true;
    options.show = show;
    geodesic_draping::plotDrapeComparisonResult(mesh, completeFieldResult, fastResult, options);

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

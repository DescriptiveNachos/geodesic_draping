#include "fixture_io.h"
#include "geodesic_draping/custom_signed_heat.h"
#include "geodesic_draping/generator_tracing.h"
#include "geodesic_draping/geometrycentral_adapter.h"
#include "geodesic_draping/seed_projection.h"

#include <algorithm>
#include <chrono>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

using Clock = std::chrono::steady_clock;

struct Summary {
  double mean = 0.0;
  double median = 0.0;
  double min = 0.0;
  double max = 0.0;
};

struct HeatRun {
  double constructSeconds = 0.0;
  double traceSeconds = 0.0;
  std::array<geodesic_draping::CustomSignedHeatStageTimings, 2> families;
};

double secondsSince(const Clock::time_point& start) {
  return std::chrono::duration<double>(Clock::now() - start).count();
}

Summary summarize(std::vector<double> values) {
  if (values.empty()) {
    throw std::runtime_error("cannot summarize an empty timing series");
  }
  std::sort(values.begin(), values.end());
  const double sum = std::accumulate(values.begin(), values.end(), 0.0);
  Summary summary;
  summary.mean = sum / static_cast<double>(values.size());
  summary.median = values[values.size() / 2];
  if (values.size() % 2 == 0) {
    summary.median = 0.5 * (values[values.size() / 2 - 1] + values[values.size() / 2]);
  }
  summary.min = values.front();
  summary.max = values.back();
  return summary;
}

geodesic_draping::SignedHeatSolveOptions fixtureHeatOptions() {
  geodesic_draping::SignedHeatSolveOptions options;
  options.preserveSourceNormals = false;
  options.levelSetConstraint = geometrycentral::LevelSetConstraint::None;
  options.softLevelSetWeight = -1.0;
  return options;
}

std::vector<double> collect(const std::vector<HeatRun>& runs, double HeatRun::*field) {
  std::vector<double> values;
  values.reserve(runs.size());
  for (const HeatRun& run : runs) {
    values.push_back(run.*field);
  }
  return values;
}

std::vector<double> collectFamily(const std::vector<HeatRun>& runs,
                                  size_t family,
                                  double geodesic_draping::CustomSignedHeatStageTimings::*field) {
  std::vector<double> values;
  values.reserve(runs.size());
  for (const HeatRun& run : runs) {
    values.push_back(run.families[family].*field);
  }
  return values;
}

std::vector<double> collectFamilySum(const std::vector<HeatRun>& runs,
                                     double geodesic_draping::CustomSignedHeatStageTimings::*field) {
  std::vector<double> values;
  values.reserve(runs.size());
  for (const HeatRun& run : runs) {
    values.push_back(run.families[0].*field + run.families[1].*field);
  }
  return values;
}

void printSummary(const std::string& scope, const std::string& stage, const Summary& summary) {
  std::cout << std::left << std::setw(12) << scope
            << std::setw(18) << stage
            << std::right << std::fixed << std::setprecision(6)
            << std::setw(12) << summary.mean
            << std::setw(12) << summary.median
            << std::setw(12) << summary.min
            << std::setw(12) << summary.max
            << "\n";
}

void printTimedHeatSummary(const std::string& label, const std::vector<HeatRun>& runs) {
  std::cout << "\n" << label << "\n";
  std::cout << std::left << std::setw(12) << "scope"
            << std::setw(18) << "stage"
            << std::right
            << std::setw(12) << "mean"
            << std::setw(12) << "median"
            << std::setw(12) << "min"
            << std::setw(12) << "max"
            << "\n";
  printSummary("setup", "construct", summarize(collect(runs, &HeatRun::constructSeconds)));
  printSummary("setup", "trace", summarize(collect(runs, &HeatRun::traceSeconds)));
  for (size_t family = 0; family < 2; ++family) {
    const std::string scope = "family_" + std::to_string(family);
    printSummary(scope, "curve_convert", summarize(collectFamily(runs, family, &geodesic_draping::CustomSignedHeatStageTimings::curveConversionSeconds)));
    printSummary(scope, "preprocess", summarize(collectFamily(runs, family, &geodesic_draping::CustomSignedHeatStageTimings::preprocessSeconds)));
    printSummary(scope, "source", summarize(collectFamily(runs, family, &geodesic_draping::CustomSignedHeatStageTimings::sourceSeconds)));
    printSummary(scope, "diffuse", summarize(collectFamily(runs, family, &geodesic_draping::CustomSignedHeatStageTimings::diffuseSeconds)));
    printSummary(scope, "normalize", summarize(collectFamily(runs, family, &geodesic_draping::CustomSignedHeatStageTimings::normalizeSeconds)));
    printSummary(scope, "distance", summarize(collectFamily(runs, family, &geodesic_draping::CustomSignedHeatStageTimings::distanceSeconds)));
    printSummary(scope, "total", summarize(collectFamily(runs, family, &geodesic_draping::CustomSignedHeatStageTimings::totalSeconds)));
  }
  printSummary("both", "curve_convert", summarize(collectFamilySum(runs, &geodesic_draping::CustomSignedHeatStageTimings::curveConversionSeconds)));
  printSummary("both", "preprocess", summarize(collectFamilySum(runs, &geodesic_draping::CustomSignedHeatStageTimings::preprocessSeconds)));
  printSummary("both", "source", summarize(collectFamilySum(runs, &geodesic_draping::CustomSignedHeatStageTimings::sourceSeconds)));
  printSummary("both", "diffuse", summarize(collectFamilySum(runs, &geodesic_draping::CustomSignedHeatStageTimings::diffuseSeconds)));
  printSummary("both", "normalize", summarize(collectFamilySum(runs, &geodesic_draping::CustomSignedHeatStageTimings::normalizeSeconds)));
  printSummary("both", "distance", summarize(collectFamilySum(runs, &geodesic_draping::CustomSignedHeatStageTimings::distanceSeconds)));
  printSummary("both", "total", summarize(collectFamilySum(runs, &geodesic_draping::CustomSignedHeatStageTimings::totalSeconds)));
}

HeatRun runColdOnce(const geodesic_draping::SurfaceMeshData& mesh,
                    const geodesic_draping::Vec2& seedXY,
                    double angleDegrees,
                    bool computeDistance) {
  HeatRun run;
  auto start = Clock::now();
  auto surface = geodesic_draping::makeGeometryCentralSurface(mesh);
  geodesic_draping::CustomSignedHeatSolver heatSolver(
      surface,
      fixtureHeatOptions().diffusionTimeCoefficient);
  run.constructSeconds = secondsSince(start);

  start = Clock::now();
  const auto seed = geodesic_draping::projectPointXYToMesh(mesh, seedXY);
  if (!seed) {
    throw std::runtime_error("seed projection failed");
  }
  const auto directions = geodesic_draping::generateOrthogonalDirections(angleDegrees);
  const auto generators = geodesic_draping::traceGenerators(surface, seed->surfacePoint, directions);
  const auto sourceCurves = geodesic_draping::pairOppositeGeneratorTraces(generators);
  run.traceSeconds = secondsSince(start);

  const auto timed = heatSolver.solveTimed(sourceCurves, fixtureHeatOptions(), computeDistance);
  run.families = {timed[0].timings, timed[1].timings};
  return run;
}

std::vector<HeatRun> benchmarkCold(const geodesic_draping::SurfaceMeshData& mesh,
                                   const geodesic_draping::Vec2& seedXY,
                                   double angleDegrees,
                                   bool computeDistance,
                                   int warmupRuns,
                                   int measuredRuns) {
  for (int i = 0; i < warmupRuns; ++i) {
    (void)runColdOnce(mesh, seedXY, angleDegrees, computeDistance);
  }
  std::vector<HeatRun> runs;
  runs.reserve(static_cast<size_t>(measuredRuns));
  for (int i = 0; i < measuredRuns; ++i) {
    runs.push_back(runColdOnce(mesh, seedXY, angleDegrees, computeDistance));
  }
  return runs;
}

std::vector<HeatRun> benchmarkWarm(const geodesic_draping::SurfaceMeshData& mesh,
                                   const geodesic_draping::Vec2& seedXY,
                                   double angleDegrees,
                                   bool computeDistance,
                                   int warmupRuns,
                                   int measuredRuns) {
  HeatRun run;
  auto start = Clock::now();
  auto surface = geodesic_draping::makeGeometryCentralSurface(mesh);
  geodesic_draping::CustomSignedHeatSolver heatSolver(
      surface,
      fixtureHeatOptions().diffusionTimeCoefficient);
  run.constructSeconds = secondsSince(start);

  start = Clock::now();
  const auto seed = geodesic_draping::projectPointXYToMesh(mesh, seedXY);
  if (!seed) {
    throw std::runtime_error("seed projection failed");
  }
  const auto directions = geodesic_draping::generateOrthogonalDirections(angleDegrees);
  const auto generators = geodesic_draping::traceGenerators(surface, seed->surfacePoint, directions);
  const auto sourceCurves = geodesic_draping::pairOppositeGeneratorTraces(generators);
  run.traceSeconds = secondsSince(start);

  for (int i = 0; i < warmupRuns; ++i) {
    (void)heatSolver.solveTimed(sourceCurves, fixtureHeatOptions(), computeDistance);
  }
  std::vector<HeatRun> runs;
  runs.reserve(static_cast<size_t>(measuredRuns));
  for (int i = 0; i < measuredRuns; ++i) {
    const auto timed = heatSolver.solveTimed(sourceCurves, fixtureHeatOptions(), computeDistance);
    HeatRun measuredRun;
    measuredRun.constructSeconds = 0.0;
    measuredRun.traceSeconds = 0.0;
    measuredRun.families = {timed[0].timings, timed[1].timings};
    runs.push_back(measuredRun);
  }
  return runs;
}

void printUsage(const char* argv0) {
  std::cout << "Usage:\n"
            << "  " << argv0 << " [fixture-name] [--warmup-runs n] [--measured-runs n]\n\n"
            << "Default fixture: demo_part\n";
}

} // namespace

int main(int argc, char** argv) {
  try {
    std::string fixtureName = "demo_part";
    int warmupRuns = 2;
    int measuredRuns = 10;

    for (int i = 1; i < argc; ++i) {
      const std::string arg = argv[i];
      if (arg == "--help" || arg == "-h") {
        printUsage(argv[0]);
        return 0;
      }
      if (arg == "--warmup-runs") {
        if (i + 1 >= argc) throw std::runtime_error("--warmup-runs requires a value");
        warmupRuns = std::stoi(argv[++i]);
        continue;
      }
      if (arg == "--measured-runs") {
        if (i + 1 >= argc) throw std::runtime_error("--measured-runs requires a value");
        measuredRuns = std::stoi(argv[++i]);
        continue;
      }
      fixtureName = arg;
    }

    const std::filesystem::path fixtureDir =
        std::filesystem::path(GEODESIC_DRAPING_TEST_DATA_DIR) / fixtureName;
    const auto mesh = geodesic_draping::fixture_io::loadMesh(fixtureDir);
    const auto seedXY = geodesic_draping::fixture_io::loadSeedXY(fixtureDir);
    const double angleDegrees = geodesic_draping::fixture_io::loadAngleDegrees(fixtureDir);

    std::cout << "fixture=" << fixtureName
              << " warmup_runs=" << warmupRuns
              << " measured_runs=" << measuredRuns << "\n";
    printTimedHeatSummary("cold fast heat stages", benchmarkCold(mesh, seedXY, angleDegrees, false, warmupRuns, measuredRuns));
    printTimedHeatSummary("cold complete heat stages", benchmarkCold(mesh, seedXY, angleDegrees, true, warmupRuns, measuredRuns));
    printTimedHeatSummary("warm fast heat stages", benchmarkWarm(mesh, seedXY, angleDegrees, false, warmupRuns, measuredRuns));
    printTimedHeatSummary("warm complete heat stages", benchmarkWarm(mesh, seedXY, angleDegrees, true, warmupRuns, measuredRuns));
    return 0;
  } catch (const std::exception& e) {
    std::cerr << "error: " << e.what() << "\n\n";
    printUsage(argv[0]);
    return 1;
  }
}

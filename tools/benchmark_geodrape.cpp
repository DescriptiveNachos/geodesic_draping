#include "fixture_io.h"
#include "geodesic_draping/field_processing.h"
#include "geodesic_draping/geodrape.h"
#include "geodesic_draping/geometrycentral_adapter.h"
#include "geodesic_draping/seed_projection.h"
#include "geodesic_draping/signed_heat.h"
#include "geodesic_draping/signed_vector_heat.h"

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

struct StepTimings {
  double constructSeconds = 0.0;
  double originSeconds = 0.0;
  double generatorsSeconds = 0.0;
  double heatSeconds = 0.0;
  double fieldsSeconds = 0.0;
  double solveSeconds = 0.0;
  double totalSeconds = 0.0;
};

struct Summary {
  double mean = 0.0;
  double median = 0.0;
  double min = 0.0;
  double max = 0.0;
};

struct BenchmarkResult {
  std::vector<StepTimings> runs;
  Summary total;
  Summary solve;
  Summary construct;
  Summary origin;
  Summary generators;
  Summary heat;
  Summary fields;
};

struct ReferenceSummary {
  double totalMean = 0.0;
  double solveMean = 0.0;
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

std::vector<double> collect(const std::vector<StepTimings>& runs, double StepTimings::*field) {
  std::vector<double> values;
  values.reserve(runs.size());
  for (const StepTimings& run : runs) {
    values.push_back(run.*field);
  }
  return values;
}

BenchmarkResult summarizeRuns(std::vector<StepTimings> runs) {
  BenchmarkResult result;
  result.runs = std::move(runs);
  result.total = summarize(collect(result.runs, &StepTimings::totalSeconds));
  result.solve = summarize(collect(result.runs, &StepTimings::solveSeconds));
  result.construct = summarize(collect(result.runs, &StepTimings::constructSeconds));
  result.origin = summarize(collect(result.runs, &StepTimings::originSeconds));
  result.generators = summarize(collect(result.runs, &StepTimings::generatorsSeconds));
  result.heat = summarize(collect(result.runs, &StepTimings::heatSeconds));
  result.fields = summarize(collect(result.runs, &StepTimings::fieldsSeconds));
  return result;
}

geodesic_draping::SignedHeatSolveOptions fixtureHeatOptions() {
  geodesic_draping::SignedHeatSolveOptions options;
  options.preserveSourceNormals = false;
  options.levelSetConstraint = geometrycentral::LevelSetConstraint::None;
  options.softLevelSetWeight = -1.0;
  return options;
}

StepTimings runCompleteOnce(const geodesic_draping::SurfaceMeshData& mesh,
                            const geodesic_draping::Vec2& seedXY,
                            double angleDegrees) {
  StepTimings timings;
  const auto constructStart = Clock::now();
  auto surface = geodesic_draping::makeGeometryCentralSurface(mesh);
  geodesic_draping::SignedHeatDistanceSolver distanceSolver(
      surface,
      fixtureHeatOptions().diffusionTimeCoefficient);
  timings.constructSeconds = secondsSince(constructStart);

  const auto solveStart = Clock::now();
  auto phaseStart = Clock::now();
  const auto seed = geodesic_draping::projectPointXYToMesh(mesh, seedXY);
  if (!seed) {
    throw std::runtime_error("seed projection failed");
  }
  const auto directions = geodesic_draping::generateOrthogonalDirections(angleDegrees);
  timings.originSeconds = secondsSince(phaseStart);

  phaseStart = Clock::now();
  const auto generators = geodesic_draping::traceGenerators(surface, seed->surfacePoint, directions);
  const auto sourceCurves = geodesic_draping::pairOppositeGeneratorTraces(generators);
  timings.generatorsSeconds = secondsSince(phaseStart);

  phaseStart = Clock::now();
  const auto distances = distanceSolver.computeDistances(sourceCurves, fixtureHeatOptions());
  timings.heatSeconds = secondsSince(phaseStart);

  phaseStart = Clock::now();
  const auto grad0 = geodesic_draping::computeVertexScalarGradients(mesh, distances[0]);
  const auto grad1 = geodesic_draping::computeVertexScalarGradients(mesh, distances[1]);
  const auto shear = geodesic_draping::computeShearAnglesDegrees(grad0, grad1);
  (void)shear;
  timings.fieldsSeconds = secondsSince(phaseStart);

  timings.solveSeconds = secondsSince(solveStart);
  timings.totalSeconds = timings.constructSeconds + timings.solveSeconds;
  return timings;
}

StepTimings runFastOnce(const geodesic_draping::SurfaceMeshData& mesh,
                        const geodesic_draping::Vec2& seedXY,
                        double angleDegrees) {
  StepTimings timings;
  const auto constructStart = Clock::now();
  auto surface = geodesic_draping::makeGeometryCentralSurface(mesh);
  geodesic_draping::CustomSignedHeatSolver customHeatSolver(
      surface,
      fixtureHeatOptions().diffusionTimeCoefficient);
  timings.constructSeconds = secondsSince(constructStart);

  const auto solveStart = Clock::now();
  auto phaseStart = Clock::now();
  const auto seed = geodesic_draping::projectPointXYToMesh(mesh, seedXY);
  if (!seed) {
    throw std::runtime_error("seed projection failed");
  }
  const auto directions = geodesic_draping::generateOrthogonalDirections(angleDegrees);
  timings.originSeconds = secondsSince(phaseStart);

  phaseStart = Clock::now();
  const auto generators = geodesic_draping::traceGenerators(surface, seed->surfacePoint, directions);
  const auto sourceCurves = geodesic_draping::pairOppositeGeneratorTraces(generators);
  timings.generatorsSeconds = secondsSince(phaseStart);

  phaseStart = Clock::now();
  std::array<geodesic_draping::CustomSignedHeatResult, 2> heats;
  for (size_t i = 0; i < heats.size(); ++i) {
    heats[i].diffusion = customHeatSolver.solveDiffusedEdgeHeatField(sourceCurves.curves[i], fixtureHeatOptions());
    heats[i].normalizedFaceDirections =
        geodesic_draping::sampleAndNormalizeFaceDirections(surface, heats[i].diffusion.diffusedEdgeHeatField);
    heats[i].vertexDirections =
        geodesic_draping::averageFaceDirectionsToVerticesReference(surface, heats[i].normalizedFaceDirections);
  }
  timings.heatSeconds = secondsSince(phaseStart);

  phaseStart = Clock::now();
  const auto shear = geodesic_draping::computeShearAnglesDegrees(heats[0].vertexDirections, heats[1].vertexDirections);
  (void)shear;
  timings.fieldsSeconds = secondsSince(phaseStart);

  timings.solveSeconds = secondsSince(solveStart);
  timings.totalSeconds = timings.constructSeconds + timings.solveSeconds;
  return timings;
}

BenchmarkResult benchmark(const geodesic_draping::SurfaceMeshData& mesh,
                          const geodesic_draping::Vec2& seedXY,
                          double angleDegrees,
                          bool fast,
                          int warmupRuns,
                          int measuredRuns) {
  for (int i = 0; i < warmupRuns; ++i) {
    if (fast) {
      runFastOnce(mesh, seedXY, angleDegrees);
    } else {
      runCompleteOnce(mesh, seedXY, angleDegrees);
    }
  }

  std::vector<StepTimings> runs;
  runs.reserve(static_cast<size_t>(measuredRuns));
  for (int i = 0; i < measuredRuns; ++i) {
    runs.push_back(fast ? runFastOnce(mesh, seedXY, angleDegrees)
                        : runCompleteOnce(mesh, seedXY, angleDegrees));
  }
  return summarizeRuns(std::move(runs));
}

double numberAfterKey(const std::string& text, const std::string& key, size_t startPos) {
  const std::string quotedKey = "\"" + key + "\"";
  const size_t keyPos = text.find(quotedKey, startPos);
  if (keyPos == std::string::npos) {
    throw std::runtime_error("missing JSON key " + key);
  }
  const size_t colon = text.find(':', keyPos + quotedKey.size());
  if (colon == std::string::npos) {
    throw std::runtime_error("missing JSON value for " + key);
  }
  const char* cursor = text.c_str() + colon + 1;
  char* end = nullptr;
  const double value = std::strtod(cursor, &end);
  if (end == cursor) {
    throw std::runtime_error("expected numeric JSON value for " + key);
  }
  return value;
}

ReferenceSummary loadReferenceSummary(const std::filesystem::path& fixtureDir, const std::string& section) {
  const std::string text = geodesic_draping::fixture_io::readText(fixtureDir / "timings_reference.json");
  const std::string quotedSection = "\"" + section + "\"";
  const size_t sectionPos = text.find(quotedSection);
  if (sectionPos == std::string::npos) {
    throw std::runtime_error("missing reference timing section " + section);
  }
  ReferenceSummary summary;
  summary.totalMean = numberAfterKey(text, "mean_seconds", sectionPos);
  summary.solveMean = numberAfterKey(text, "solve_mean_seconds", sectionPos);
  return summary;
}

void printPath(const std::string& fixtureName,
               const std::string& pathName,
               const BenchmarkResult& cpp,
               const ReferenceSummary& ref) {
  const double totalSpeedup = ref.totalMean / cpp.total.mean;
  const double solveSpeedup = ref.solveMean / cpp.solve.mean;
  std::cout << std::left << std::setw(13) << fixtureName
            << std::setw(10) << pathName
            << std::right << std::fixed << std::setprecision(6)
            << std::setw(12) << cpp.total.mean
            << std::setw(12) << ref.totalMean
            << std::setw(10) << totalSpeedup
            << std::setw(12) << cpp.solve.mean
            << std::setw(12) << ref.solveMean
            << std::setw(10) << solveSpeedup
            << std::setw(12) << cpp.construct.mean
            << std::setw(12) << cpp.origin.mean
            << std::setw(12) << cpp.generators.mean
            << std::setw(12) << cpp.heat.mean
            << std::setw(12) << cpp.fields.mean
            << "\n";
}

void printUsage(const char* argv0) {
  std::cout << "Usage:\n"
            << "  " << argv0 << " [fixture-name|all] [--warmup-runs n] [--measured-runs n]\n\n"
            << "Default fixture: all\n";
}

} // namespace

int main(int argc, char** argv) {
  try {
    std::vector<std::string> fixtures = {"tiny_planar", "small_curved", "demo_part"};
    int warmupRuns = 1;
    int measuredRuns = 5;

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
      if (arg != "all") {
        fixtures = {arg};
      }
    }

    if (warmupRuns < 0 || measuredRuns <= 0) {
      throw std::runtime_error("warmup runs must be non-negative and measured runs must be positive");
    }

    const std::filesystem::path root = GEODESIC_DRAPING_TEST_DATA_DIR;
    std::cout << "warmup_runs=" << warmupRuns << " measured_runs=" << measuredRuns << "\n";
    std::cout << std::left << std::setw(13) << "fixture"
              << std::setw(10) << "path"
              << std::right
              << std::setw(12) << "cpp_total"
              << std::setw(12) << "ref_total"
              << std::setw(10) << "x_total"
              << std::setw(12) << "cpp_solve"
              << std::setw(12) << "ref_solve"
              << std::setw(10) << "x_solve"
              << std::setw(12) << "construct"
              << std::setw(12) << "origin"
              << std::setw(12) << "generators"
              << std::setw(12) << "heat"
              << std::setw(12) << "fields"
              << "\n";

    for (const std::string& fixtureName : fixtures) {
      const std::filesystem::path fixtureDir = root / fixtureName;
      const auto mesh = geodesic_draping::fixture_io::loadMesh(fixtureDir);
      const auto seedXY = geodesic_draping::fixture_io::loadSeedXY(fixtureDir);
      const double angleDegrees = geodesic_draping::fixture_io::loadAngleDegrees(fixtureDir);

      const BenchmarkResult complete = benchmark(mesh, seedXY, angleDegrees, false, warmupRuns, measuredRuns);
      const BenchmarkResult fast = benchmark(mesh, seedXY, angleDegrees, true, warmupRuns, measuredRuns);
      printPath(fixtureName, "complete", complete, loadReferenceSummary(fixtureDir, "complete"));
      printPath(fixtureName, "fast", fast, loadReferenceSummary(fixtureDir, "fast"));
    }
    return 0;
  } catch (const std::exception& e) {
    std::cerr << "error: " << e.what() << "\n\n";
    printUsage(argv[0]);
    return 1;
  }
}

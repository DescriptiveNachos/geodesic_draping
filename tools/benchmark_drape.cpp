#include "fixture_io.h"
#include "geodesic_draping/geodrape.h"
#include "geodesic_draping/strings.h"

#include <algorithm>
#include <chrono>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

using Clock = std::chrono::steady_clock;

volatile size_t gResultSink = 0;

struct Options {
  std::string fixtureName = "demo_part";
  size_t warmupRuns = 1;
  size_t measuredRuns = 10;
  geodesic_draping::RetrievalDomain retrievalDomain =
      geodesic_draping::RetrievalDomain::Extrinsic;
};

struct TimingStats {
  double meanMs = 0.0;
  double medianMs = 0.0;
  double minMs = 0.0;
  double maxMs = 0.0;
};

void printUsage(const char* argv0) {
  std::cout << "Usage:\n"
            << "  " << argv0 << " [fixture-name]\n"
            << "    [--warmup-runs n]\n"
            << "    [--measured-runs n]\n"
            << "    [--domain extrinsic|subdivision]\n\n"
            << "Default fixture: demo_part\n"
            << "Default domain: extrinsic\n";
}

geodesic_draping::RetrievalDomain parseDomain(const std::string& value) {
  if (value == "extrinsic" || value == "subdivision") {
    return geodesic_draping::parseRetrievalDomain(value);
  }
  throw std::runtime_error("--domain must be one of: extrinsic, subdivision");
}

Options parseArgs(int argc, char** argv) {
  Options options;
  for (int i = 1; i < argc; ++i) {
    const std::string arg = argv[i];
    if (arg == "--help" || arg == "-h") {
      printUsage(argv[0]);
      std::exit(0);
    }
    if (arg == "--warmup-runs") {
      if (++i >= argc) throw std::runtime_error("--warmup-runs requires a value");
      options.warmupRuns = static_cast<size_t>(std::stoull(argv[i]));
    } else if (arg == "--measured-runs") {
      if (++i >= argc) throw std::runtime_error("--measured-runs requires a value");
      options.measuredRuns = static_cast<size_t>(std::stoull(argv[i]));
    } else if (arg == "--domain") {
      if (++i >= argc) throw std::runtime_error("--domain requires a value");
      options.retrievalDomain = parseDomain(argv[i]);
    } else if (arg.rfind("--", 0) == 0) {
      throw std::runtime_error("unknown option: " + arg);
    } else {
      options.fixtureName = arg;
    }
  }
  if (options.measuredRuns == 0) {
    throw std::runtime_error("--measured-runs must be greater than zero");
  }
  return options;
}

void consumeResult(const geodesic_draping::DrapeResult& result) {
  size_t value = result.mesh ? result.mesh->nVertices() + result.mesh->nFaces() : 0;
  if (result.distances) value += 17;
  if (result.faceShear) value += 31;
  if (result.vertexShear) value += 43;
  if (result.directionFields) value += 59;
  gResultSink += value;
}

TimingStats stats(std::vector<double> samplesMs) {
  std::sort(samplesMs.begin(), samplesMs.end());
  TimingStats s;
  s.minMs = samplesMs.front();
  s.maxMs = samplesMs.back();
  s.medianMs = samplesMs[samplesMs.size() / 2];
  double sum = 0.0;
  for (double value : samplesMs) {
    sum += value;
  }
  s.meanMs = sum / static_cast<double>(samplesMs.size());
  return s;
}

template <typename Fn>
std::vector<double> measure(size_t runs, Fn&& fn) {
  std::vector<double> samples;
  samples.reserve(runs);
  for (size_t i = 0; i < runs; ++i) {
    const Clock::time_point start = Clock::now();
    fn();
    const Clock::time_point stop = Clock::now();
    samples.push_back(
        std::chrono::duration<double, std::milli>(stop - start).count());
  }
  return samples;
}

std::vector<double> measureCold(
    const geodesic_draping::SurfaceMeshData& meshData,
    const geodesic_draping::Vec2& seedXY,
    double fabricAngle,
    const geodesic_draping::SignedHeatSolveOptions& heatOptions,
    const geodesic_draping::DrapeSolveOptions& solveOptions,
    const geodesic_draping::RetrievalOptions& retrievalOptions,
    size_t warmupRuns,
    size_t measuredRuns) {
  for (size_t i = 0; i < warmupRuns; ++i) {
    consumeResult(geodesic_draping::solveDrape(
        meshData, seedXY, fabricAngle, heatOptions, solveOptions, retrievalOptions));
  }
  return measure(measuredRuns, [&]() {
    consumeResult(geodesic_draping::solveDrape(
        meshData, seedXY, fabricAngle, heatOptions, solveOptions, retrievalOptions));
  });
}

std::vector<double> measureWarm(
    const geodesic_draping::SurfaceMeshData& meshData,
    const geodesic_draping::Vec2& seedXY,
    double fabricAngle,
    const geodesic_draping::SignedHeatSolveOptions& heatOptions,
    const geodesic_draping::DrapeSolveOptions& solveOptions,
    const geodesic_draping::RetrievalOptions& retrievalOptions,
    size_t warmupRuns,
    size_t measuredRuns) {
  geodesic_draping::GeoDrapeSolver solver(meshData, heatOptions);
  for (size_t i = 0; i < warmupRuns; ++i) {
    consumeResult(solver.solve(seedXY, fabricAngle, solveOptions, retrievalOptions));
  }
  return measure(measuredRuns, [&]() {
    consumeResult(solver.solve(seedXY, fabricAngle, solveOptions, retrievalOptions));
  });
}

void printRow(const std::string& mode,
              const std::string& scope,
              const TimingStats& stats) {
  std::cout << std::left << std::setw(10) << mode
            << std::setw(8) << scope
            << std::right << std::fixed << std::setprecision(3)
            << std::setw(12) << stats.meanMs
            << std::setw(12) << stats.medianMs
            << std::setw(12) << stats.minMs
            << std::setw(12) << stats.maxMs << "\n";
}

} // namespace

int main(int argc, char** argv) {
  try {
    const Options options = parseArgs(argc, argv);
    const std::filesystem::path fixtureDir =
        std::filesystem::path(GEODESIC_DRAPING_TEST_DATA_DIR) / options.fixtureName;
    const geodesic_draping::SurfaceMeshData meshData =
        geodesic_draping::fixture_io::loadMesh(fixtureDir);
    const geodesic_draping::Vec2 seedXY =
        geodesic_draping::fixture_io::loadSeedXY(fixtureDir);
    const double fabricAngle =
        geodesic_draping::fixture_io::loadAngleDegrees(fixtureDir);
    const geodesic_draping::SignedHeatSolveOptions heatOptions;

    geodesic_draping::RetrievalOptions retrievalOptions;
    retrievalOptions.domain = options.retrievalDomain;

    std::cout << "fixture: " << options.fixtureName
              << "  vertices: " << meshData.vertices.size()
              << "  faces: " << meshData.faces.size()
              << "  domain: " << geodesic_draping::retrievalDomainName(options.retrievalDomain)
              << "  warmups: " << options.warmupRuns
              << "  measured: " << options.measuredRuns << "\n\n";
    std::cout << std::left << std::setw(10) << "mode"
              << std::setw(8) << "scope"
              << std::right
              << std::setw(12) << "mean_ms"
              << std::setw(12) << "median_ms"
              << std::setw(12) << "min_ms"
              << std::setw(12) << "max_ms" << "\n";

    const std::vector<geodesic_draping::DrapeSolveMode> modes{
        geodesic_draping::DrapeSolveMode::Fast,
        geodesic_draping::DrapeSolveMode::Hybrid,
        geodesic_draping::DrapeSolveMode::Complete,
    };

    for (geodesic_draping::DrapeSolveMode mode : modes) {
      geodesic_draping::DrapeSolveOptions solveOptions;
      solveOptions.mode = mode;
      const std::string name = geodesic_draping::drapeSolveModeName(mode);
      printRow(name, "cold", stats(measureCold(
                         meshData,
                         seedXY,
                         fabricAngle,
                         heatOptions,
                         solveOptions,
                         retrievalOptions,
                         options.warmupRuns,
                         options.measuredRuns)));
      printRow(name, "warm", stats(measureWarm(
                         meshData,
                         seedXY,
                         fabricAngle,
                         heatOptions,
                         solveOptions,
                         retrievalOptions,
                         options.warmupRuns,
                         options.measuredRuns)));
    }
  } catch (const std::exception& e) {
    std::cerr << "benchmark_drape failed: " << e.what() << "\n\n";
    printUsage(argv[0]);
    return 1;
  }
  return 0;
}

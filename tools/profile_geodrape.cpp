#include "fixture_io.h"
#include "geodesic_draping/custom_signed_heat.h"
#include "geodesic_draping/field_processing.h"
#include "geodesic_draping/geodrape.h"
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

std::string modeName(geodesic_draping::DrapeSolveMode mode) {
  switch (mode) {
  case geodesic_draping::DrapeSolveMode::Fast:
    return "fast";
  case geodesic_draping::DrapeSolveMode::Hybrid:
    return "hybrid";
  case geodesic_draping::DrapeSolveMode::Complete:
    return "complete";
  }
  return "unknown";
}

std::string domainName(geodesic_draping::ResultDomain domain) {
  switch (domain) {
  case geodesic_draping::ResultDomain::Intrinsic:
    return "intrinsic";
  case geodesic_draping::ResultDomain::Extrinsic:
    return "extrinsic";
  case geodesic_draping::ResultDomain::Subdivision:
    return "subdivision";
  }
  return "unknown";
}

std::string backendName(geodesic_draping::IntrinsicTriangulationBackend backend) {
  switch (backend) {
  case geodesic_draping::IntrinsicTriangulationBackend::Signpost:
    return "signpost";
  case geodesic_draping::IntrinsicTriangulationBackend::IntegerCoordinates:
    return "integer";
  }
  return "unknown";
}

void consume(const geodesic_draping::DrapeResult& result) {
  volatile size_t sink = result.mesh.faces.size();
  if (result.faceShear) {
    sink += result.faceShear->size();
  }
  if (result.vertexShear) {
    sink += result.vertexShear->size();
  }
  if (result.distances) {
    sink += (*result.distances)[0].size() + (*result.distances)[1].size();
  }
  (void)sink;
}

geodesic_draping::DrapeSolveOptions solveOptions(geodesic_draping::DrapeSolveMode mode,
                                                 geodesic_draping::ResultDomain retrieval,
                                                 bool sampleVertexShear) {
  geodesic_draping::DrapeSolveOptions options;
  options.mode = mode;
  options.retrieval = retrieval;
  options.sampleVertexShear = sampleVertexShear;
  return options;
}

Summary timeSolverConstruction(const geodesic_draping::SurfaceMeshData& mesh,
                               const geodesic_draping::IntrinsicConstructionOptions& intrinsic,
                               const geodesic_draping::RefinementOptions& refinement,
                               int warmupRuns,
                               int measuredRuns) {
  for (int i = 0; i < warmupRuns; ++i) {
    geodesic_draping::GeoDrapeSolver solver(mesh, fixtureHeatOptions(), intrinsic, refinement);
    (void)solver;
  }

  std::vector<double> values;
  values.reserve(static_cast<size_t>(measuredRuns));
  for (int i = 0; i < measuredRuns; ++i) {
    const auto start = Clock::now();
    geodesic_draping::GeoDrapeSolver solver(mesh, fixtureHeatOptions(), intrinsic, refinement);
    (void)solver;
    values.push_back(secondsSince(start));
  }
  return summarize(std::move(values));
}

Summary timeSolve(const geodesic_draping::SurfaceMeshData& mesh,
                  const geodesic_draping::Vec2& seedXY,
                  double angleDegrees,
                  const geodesic_draping::IntrinsicConstructionOptions& intrinsic,
                  const geodesic_draping::RefinementOptions& refinement,
                  const geodesic_draping::DrapeSolveOptions& options,
                  int warmupRuns,
                  int measuredRuns) {
  geodesic_draping::GeoDrapeSolver solver(mesh, fixtureHeatOptions(), intrinsic, refinement);
  for (int i = 0; i < warmupRuns; ++i) {
    consume(solver.solve(seedXY, angleDegrees, options));
  }

  std::vector<double> values;
  values.reserve(static_cast<size_t>(measuredRuns));
  for (int i = 0; i < measuredRuns; ++i) {
    const auto start = Clock::now();
    consume(solver.solve(seedXY, angleDegrees, options));
    values.push_back(secondsSince(start));
  }
  return summarize(std::move(values));
}

Summary timeRetrieve(const geodesic_draping::SurfaceMeshData& mesh,
                     const geodesic_draping::Vec2& seedXY,
                     double angleDegrees,
                     const geodesic_draping::IntrinsicConstructionOptions& intrinsic,
                     const geodesic_draping::RefinementOptions& refinement,
                     geodesic_draping::DrapeSolveMode mode,
                     geodesic_draping::ResultDomain retrieval,
                     bool sampleVertexShear,
                     int warmupRuns,
                     int measuredRuns) {
  geodesic_draping::GeoDrapeSolver solver(mesh, fixtureHeatOptions(), intrinsic, refinement);
  consume(solver.solve(seedXY, angleDegrees, solveOptions(mode, geodesic_draping::ResultDomain::Intrinsic, false)));
  for (int i = 0; i < warmupRuns; ++i) {
    consume(solver.retrieve(retrieval, sampleVertexShear));
  }

  std::vector<double> values;
  values.reserve(static_cast<size_t>(measuredRuns));
  for (int i = 0; i < measuredRuns; ++i) {
    const auto start = Clock::now();
    consume(solver.retrieve(retrieval, sampleVertexShear));
    values.push_back(secondsSince(start));
  }
  return summarize(std::move(values));
}

Summary timeDirectExtrinsicSolve(const geodesic_draping::SurfaceMeshData& mesh,
                                 const geodesic_draping::Vec2& seedXY,
                                 double angleDegrees,
                                 geodesic_draping::DrapeSolveMode mode,
                                 int warmupRuns,
                                 int measuredRuns) {
  auto surface = geodesic_draping::makeGeometryCentralSurface(mesh);
  geodesic_draping::CustomSignedHeatSolver heatSolver(
      surface,
      fixtureHeatOptions().diffusionTimeCoefficient);

  auto runOnce = [&]() {
    const auto seed = geodesic_draping::projectPointXYToMesh(mesh, seedXY);
    if (!seed) {
      throw std::runtime_error("seed projection failed");
    }
    const auto directions = geodesic_draping::generateOrthogonalDirections(angleDegrees);
    const auto generators = geodesic_draping::traceGenerators(surface, seed->surfacePoint, directions);
    const auto sourceCurves = geodesic_draping::pairOppositeGeneratorTraces(generators);
    if (mode == geodesic_draping::DrapeSolveMode::Fast) {
      std::array<geodesic_draping::CustomSignedHeatResult, 2> heats;
      for (size_t i = 0; i < heats.size(); ++i) {
        heats[i].diffusion = heatSolver.solveDiffusedEdgeHeatField(sourceCurves.curves[i], fixtureHeatOptions());
        heats[i].normalizedFaceDirections =
            geodesic_draping::sampleAndNormalizeFaceDirections(surface, heats[i].diffusion.diffusedEdgeHeatField);
      }
      const auto faceShear = geodesic_draping::computeFaceShearAnglesDegrees(
          surface,
          heats[0].normalizedFaceDirections,
          heats[1].normalizedFaceDirections);
      consume(geodesic_draping::DrapeResult{});
      volatile size_t sink = faceShear.size();
      (void)sink;
    } else {
      const auto heatSolves = heatSolver.solve(sourceCurves, fixtureHeatOptions(), true);
      std::array<std::vector<double>, 2> distances{heatSolves[0].distance, heatSolves[1].distance};
      const auto grad0 = geodesic_draping::computeVertexScalarGradients(mesh, distances[0]);
      const auto grad1 = geodesic_draping::computeVertexScalarGradients(mesh, distances[1]);
      const auto shear = geodesic_draping::computeShearAnglesDegrees(grad0, grad1);
      volatile size_t sink = shear.size();
      (void)sink;
    }
  };

  for (int i = 0; i < warmupRuns; ++i) {
    runOnce();
  }

  std::vector<double> values;
  values.reserve(static_cast<size_t>(measuredRuns));
  for (int i = 0; i < measuredRuns; ++i) {
    const auto start = Clock::now();
    runOnce();
    values.push_back(secondsSince(start));
  }
  return summarize(std::move(values));
}

void printSummary(const std::string& kind,
                  const std::string& fixture,
                  const std::string& label,
                  const Summary& summary) {
  std::cout << std::left << std::setw(12) << kind
            << std::setw(12) << fixture
            << std::setw(24) << label
            << std::right << std::fixed << std::setprecision(6)
            << std::setw(12) << summary.mean
            << std::setw(12) << summary.median
            << std::setw(12) << summary.min
            << std::setw(12) << summary.max
            << "\n";
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
        if (i + 1 >= argc) {
          throw std::runtime_error("--warmup-runs requires a value");
        }
        warmupRuns = std::stoi(argv[++i]);
        continue;
      }
      if (arg == "--measured-runs") {
        if (i + 1 >= argc) {
          throw std::runtime_error("--measured-runs requires a value");
        }
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

    geodesic_draping::IntrinsicConstructionOptions signpost;
    geodesic_draping::IntrinsicConstructionOptions integer;
    integer.backend = geodesic_draping::IntrinsicTriangulationBackend::IntegerCoordinates;

    geodesic_draping::RefinementOptions none;
    geodesic_draping::RefinementOptions flip;
    flip.mode = geodesic_draping::RefinementMode::DelaunayFlip;
    geodesic_draping::RefinementOptions refine;
    refine.mode = geodesic_draping::RefinementMode::DelaunayRefine;
    refine.angleThreshold = 25.0;

    std::cout << "warmup_runs=" << warmupRuns << " measured_runs=" << measuredRuns << "\n";
    std::cout << std::left << std::setw(12) << "kind"
              << std::setw(12) << "fixture"
              << std::setw(24) << "label"
              << std::right
              << std::setw(12) << "mean"
              << std::setw(12) << "median"
              << std::setw(12) << "min"
              << std::setw(12) << "max"
              << "\n";

    for (const auto& config : {
             std::tuple<std::string, geodesic_draping::IntrinsicConstructionOptions, geodesic_draping::RefinementOptions>{
                 backendName(signpost.backend) + "/none", signpost, none},
             std::tuple<std::string, geodesic_draping::IntrinsicConstructionOptions, geodesic_draping::RefinementOptions>{
                 backendName(signpost.backend) + "/flip", signpost, flip},
             std::tuple<std::string, geodesic_draping::IntrinsicConstructionOptions, geodesic_draping::RefinementOptions>{
                 backendName(integer.backend) + "/none", integer, none},
             std::tuple<std::string, geodesic_draping::IntrinsicConstructionOptions, geodesic_draping::RefinementOptions>{
                 backendName(integer.backend) + "/refine", integer, refine},
         }) {
      const auto& label = std::get<0>(config);
      const auto& intrinsic = std::get<1>(config);
      const auto& refinement = std::get<2>(config);
      printSummary("construct", fixtureName, label,
                   timeSolverConstruction(mesh, intrinsic, refinement, warmupRuns, measuredRuns));
    }

    for (geodesic_draping::DrapeSolveMode mode :
         {geodesic_draping::DrapeSolveMode::Fast, geodesic_draping::DrapeSolveMode::Complete}) {
      printSummary("direct", fixtureName, modeName(mode) + "/prefactored",
                   timeDirectExtrinsicSolve(mesh, seedXY, angleDegrees, mode, warmupRuns, measuredRuns));
      for (bool sample : {false, true}) {
        const std::string sampleLabel = sample ? "+vertex" : "face-only";
        for (geodesic_draping::ResultDomain retrieval :
             {geodesic_draping::ResultDomain::Intrinsic, geodesic_draping::ResultDomain::Extrinsic}) {
          const std::string label =
              modeName(mode) + "/" + domainName(retrieval) + "/" + sampleLabel;
          printSummary("solve", fixtureName, label,
                       timeSolve(mesh, seedXY, angleDegrees, signpost, none,
                                 solveOptions(mode, retrieval, sample), warmupRuns, measuredRuns));
        }
      }
    }

    for (geodesic_draping::DrapeSolveMode mode :
         {geodesic_draping::DrapeSolveMode::Fast, geodesic_draping::DrapeSolveMode::Complete}) {
      for (bool sample : {false, true}) {
        const std::string sampleLabel = sample ? "+vertex" : "face-only";
        for (geodesic_draping::ResultDomain retrieval :
             {geodesic_draping::ResultDomain::Intrinsic, geodesic_draping::ResultDomain::Extrinsic}) {
          const std::string label =
              modeName(mode) + "/" + domainName(retrieval) + "/" + sampleLabel;
          printSummary("retrieve", fixtureName, label,
                       timeRetrieve(mesh, seedXY, angleDegrees, signpost, none,
                                    mode, retrieval, sample, warmupRuns, measuredRuns));
        }
      }
    }

    printSummary("solve", fixtureName, "complete/subdivision/flip",
                 timeSolve(mesh, seedXY, angleDegrees, signpost, flip,
                           solveOptions(geodesic_draping::DrapeSolveMode::Complete,
                                        geodesic_draping::ResultDomain::Subdivision,
                                        true),
                           warmupRuns, measuredRuns));
    printSummary("solve", fixtureName, "complete/subdivision/int-refine",
                 timeSolve(mesh, seedXY, angleDegrees, integer, refine,
                           solveOptions(geodesic_draping::DrapeSolveMode::Complete,
                                        geodesic_draping::ResultDomain::Subdivision,
                                        true),
                           warmupRuns, measuredRuns));

    return 0;
  } catch (const std::exception& e) {
    std::cerr << "error: " << e.what() << "\n\n";
    printUsage(argv[0]);
    return 1;
  }
}

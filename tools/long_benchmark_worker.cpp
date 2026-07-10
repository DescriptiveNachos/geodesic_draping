#include "fixture_io.h"
#include "geodesic_draping/geodrape.h"

#include <chrono>
#include <filesystem>
#include <iostream>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

using Clock = std::chrono::steady_clock;

double secondsSince(const Clock::time_point& start) {
  return std::chrono::duration<double>(Clock::now() - start).count();
}

geodesic_draping::SignedHeatSolveOptions fixtureHeatOptions() {
  geodesic_draping::SignedHeatSolveOptions options;
  options.preserveSourceNormals = false;
  options.levelSetConstraint = geometrycentral::LevelSetConstraint::None;
  options.softLevelSetWeight = -1.0;
  return options;
}

geodesic_draping::DrapeSolveMode parseMode(const std::string& value) {
  if (value == "fast") {
    return geodesic_draping::DrapeSolveMode::Fast;
  }
  if (value == "complete") {
    return geodesic_draping::DrapeSolveMode::Complete;
  }
  throw std::runtime_error("unknown mode: " + value);
}

geodesic_draping::DrapeSolveOptions solveOptions(const std::string& modeName) {
  geodesic_draping::DrapeSolveOptions options;
  options.mode = parseMode(modeName);
  options.retrieval = geodesic_draping::ResultDomain::Extrinsic;
  options.sampleVertexShear = true;
  return options;
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

struct ModeState {
  std::unique_ptr<geodesic_draping::GeoDrapeSolver> solver;
};

struct WorkerState {
  geodesic_draping::SurfaceMeshData mesh;
  geodesic_draping::Vec2 seedXY;
  double angleDegrees = 0.0;
  ModeState fast;
  ModeState complete;
};

ModeState& stateForMode(WorkerState& state, const std::string& modeName) {
  if (modeName == "fast") {
    return state.fast;
  }
  if (modeName == "complete") {
    return state.complete;
  }
  throw std::runtime_error("unknown mode: " + modeName);
}

void runCold(WorkerState& state, const std::string& modeName, const std::string& iteration) {
  const auto constructStart = Clock::now();
  geodesic_draping::GeoDrapeSolver solver(state.mesh, fixtureHeatOptions());
  const double constructSeconds = secondsSince(constructStart);

  const auto solveStart = Clock::now();
  consume(solver.solve(state.seedXY, state.angleDegrees, solveOptions(modeName)));
  const double solveSeconds = secondsSince(solveStart);
  const double totalSeconds = constructSeconds + solveSeconds;

  std::cout << "RESULT,cpp,cold," << modeName << "," << iteration << ","
            << constructSeconds << "," << solveSeconds << "," << totalSeconds << "\n";
}

void prepareWarm(WorkerState& state, const std::string& modeName) {
  ModeState& modeState = stateForMode(state, modeName);
  modeState.solver = std::make_unique<geodesic_draping::GeoDrapeSolver>(
      state.mesh,
      fixtureHeatOptions());
  consume(modeState.solver->solve(state.seedXY, state.angleDegrees, solveOptions(modeName)));
  std::cout << "READY," << modeName << "\n";
}

void runWarm(WorkerState& state, const std::string& modeName, const std::string& iteration) {
  ModeState& modeState = stateForMode(state, modeName);
  if (!modeState.solver) {
    prepareWarm(state, modeName);
  }

  const auto solveStart = Clock::now();
  consume(modeState.solver->solve(state.seedXY, state.angleDegrees, solveOptions(modeName)));
  const double solveSeconds = secondsSince(solveStart);

  std::cout << "RESULT,cpp,warm," << modeName << "," << iteration << ","
            << 0.0 << "," << solveSeconds << "," << solveSeconds << "\n";
}

std::vector<std::string> splitCsvLine(const std::string& line) {
  std::vector<std::string> parts;
  std::stringstream stream(line);
  std::string part;
  while (std::getline(stream, part, ',')) {
    parts.push_back(part);
  }
  return parts;
}

void printUsage(const char* argv0) {
  std::cerr << "Usage: " << argv0 << " [fixture-name]\n";
}

} // namespace

int main(int argc, char** argv) {
  try {
    const std::string fixtureName = argc > 1 ? argv[1] : "demo_part";
    const std::filesystem::path fixtureDir =
        std::filesystem::path(GEODESIC_DRAPING_TEST_DATA_DIR) / fixtureName;

    WorkerState state;
    state.mesh = geodesic_draping::fixture_io::loadMesh(fixtureDir);
    state.seedXY = geodesic_draping::fixture_io::loadSeedXY(fixtureDir);
    state.angleDegrees = geodesic_draping::fixture_io::loadAngleDegrees(fixtureDir);

    std::cout << "READY,worker\n";
    std::cout.flush();

    std::string line;
    while (std::getline(std::cin, line)) {
      if (line == "EXIT") {
        std::cout << "BYE\n";
        std::cout.flush();
        return 0;
      }

      const std::vector<std::string> parts = splitCsvLine(line);
      if (parts.size() < 3) {
        throw std::runtime_error("invalid worker command: " + line);
      }

      if (parts[0] == "PREPARE" && parts.size() == 3 && parts[1] == "warm") {
        prepareWarm(state, parts[2]);
      } else if (parts[0] == "RUN" && parts.size() == 4 && parts[1] == "cold") {
        runCold(state, parts[2], parts[3]);
      } else if (parts[0] == "RUN" && parts.size() == 4 && parts[1] == "warm") {
        runWarm(state, parts[2], parts[3]);
      } else {
        throw std::runtime_error("invalid worker command: " + line);
      }
      std::cout.flush();
    }
    return 0;
  } catch (const std::exception& e) {
    std::cerr << "error: " << e.what() << "\n";
    printUsage(argv[0]);
    return 1;
  }
}

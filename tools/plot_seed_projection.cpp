#include "fixture_io.h"
#include "geodesic_draping/plotting.h"
#include "geodesic_draping/seed_projection.h"

#include <array>
#include <cmath>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

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

void printVectorComparison(const std::string& label, const Vec3& actual, const Vec3& golden) {
  std::cout << std::left << std::setw(18) << label << " actual [" << actual.transpose() << "]"
            << "  golden [" << golden.transpose() << "]"
            << "  max_abs_diff " << maxAbsDiff(actual, golden) << "\n";
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

    const auto projection = geodesic_draping::projectPointXYToMesh(mesh, seedXY);
    if (!projection) {
      throw std::runtime_error("seed projection failed for fixture " + fixtureName);
    }

    const auto directions = geodesic_draping::generateOrthogonalDirections(angleDegrees);

    const Vec3 goldenOrigin = geodesic_draping::fixture_io::loadGoldenOrigin(fixtureDir);
    const size_t goldenFaceIndex = geodesic_draping::fixture_io::loadGoldenSeedFaceIndex(fixtureDir);
    const Vec3 goldenBarycentric = geodesic_draping::fixture_io::loadGoldenSeedBarycentric(fixtureDir);
    const std::vector<Vec3> goldenDirections = geodesic_draping::fixture_io::loadGoldenDirections(fixtureDir);

    std::cout << std::setprecision(17);
    std::cout << "Fixture: " << fixtureName << "\n"
              << "Vertices: " << mesh.vertices.size() << "  Faces: " << mesh.faces.size() << "\n"
              << "Seed XY: [" << seedXY.transpose() << "]  Angle degrees: " << angleDegrees << "\n\n";

    printVectorComparison("origin", projection->cartesian, goldenOrigin);

    std::cout << std::left << std::setw(18) << "face_index"
              << " actual " << projection->surfacePoint.faceIndex
              << "  golden " << goldenFaceIndex;
    if (projection->surfacePoint.faceIndex != goldenFaceIndex) {
      std::cout << "  note: valid tie-order differences can occur on edges/vertices";
    }
    std::cout << "\n";

    printVectorComparison("barycentric", projection->surfacePoint.barycentric, goldenBarycentric);
    std::cout << "\n";

    for (size_t i = 0; i < directions.size(); ++i) {
      printVectorComparison("direction " + std::to_string(i), directions[i], goldenDirections[i]);
    }

    ProjectionPlotOptions options;
    options.name = fixtureName + " seed projection";
    options.directionLength = directionLength;
    options.clearExisting = true;
    options.show = show;

    geodesic_draping::plotSeedProjectionStep(mesh, *projection, directions, options);

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

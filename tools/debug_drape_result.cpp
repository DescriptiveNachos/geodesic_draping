#include "fixture_io.h"
#include "geodesic_draping/geodrape.h"

#include "polyscope/curve_network.h"
#include "polyscope/point_cloud.h"
#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"

#include <array>
#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

namespace gcs = geometrycentral::surface;

using Point3 = std::array<double, 3>;
using Face3 = std::array<size_t, 3>;

void printUsage(const char* argv0) {
  std::cout << "Usage:\n"
            << "  " << argv0 << " [fixture-name]\n"
            << "    [--mode fast|hybrid|complete]\n"
            << "    [--domain extrinsic|subdivision]\n"
            << "    [--backend signpost|integer]\n"
            << "    [--refinement none|flip|refine]\n"
            << "    [--sample-vertex-shear]\n\n"
            << "Default fixture: demo_part\n"
            << "Default mode/domain: complete/subdivision\n";
}

Point3 toPoint(const geometrycentral::Vector3& p) {
  return {p.x, p.y, p.z};
}

geodesic_draping::DrapeSolveMode parseMode(const std::string& value) {
  if (value == "fast") return geodesic_draping::DrapeSolveMode::Fast;
  if (value == "hybrid") return geodesic_draping::DrapeSolveMode::Hybrid;
  if (value == "complete") return geodesic_draping::DrapeSolveMode::Complete;
  throw std::runtime_error("--mode must be one of: fast, hybrid, complete");
}

geodesic_draping::RetrievalDomain parseDomain(const std::string& value) {
  if (value == "extrinsic") return geodesic_draping::RetrievalDomain::Extrinsic;
  if (value == "subdivision") return geodesic_draping::RetrievalDomain::Subdivision;
  throw std::runtime_error("--domain must be one of: extrinsic, subdivision");
}

geodesic_draping::IntrinsicTriangulationBackend parseBackend(const std::string& value) {
  if (value == "signpost") return geodesic_draping::IntrinsicTriangulationBackend::Signpost;
  if (value == "integer" || value == "integer-coordinates") {
    return geodesic_draping::IntrinsicTriangulationBackend::IntegerCoordinates;
  }
  throw std::runtime_error("--backend must be one of: signpost, integer");
}

geodesic_draping::RefinementMode parseRefinement(const std::string& value) {
  if (value == "none") return geodesic_draping::RefinementMode::None;
  if (value == "flip") return geodesic_draping::RefinementMode::DelaunayFlip;
  if (value == "refine") return geodesic_draping::RefinementMode::DelaunayRefine;
  throw std::runtime_error("--refinement must be one of: none, flip, refine");
}

std::vector<Face3> faceIndices(const geodesic_draping::DrapeResult& result) {
  if (!result.mesh) {
    throw std::runtime_error("result does not contain a mesh");
  }
  auto& mesh = const_cast<gcs::SurfaceMesh&>(*result.mesh);
  std::vector<Face3> faces;
  faces.reserve(mesh.nFaces());
  for (gcs::Face face : mesh.faces()) {
    Face3 indices{};
    size_t i = 0;
    for (gcs::Vertex vertex : face.adjacentVertices()) {
      if (i >= indices.size()) {
        throw std::runtime_error("debug viewer expects triangular meshes");
      }
      indices[i++] = vertex.getIndex();
    }
    if (i != indices.size()) {
      throw std::runtime_error("debug viewer expects triangular meshes");
    }
    faces.push_back(indices);
  }
  return faces;
}

std::vector<Point3> vertexPositions(const geodesic_draping::DrapeResult& result) {
  if (!result.mesh) {
    throw std::runtime_error("result does not contain a mesh");
  }
  auto& mesh = const_cast<gcs::SurfaceMesh&>(*result.mesh);
  std::vector<Point3> positions;
  positions.reserve(mesh.nVertices());

  if (result.vertexPositions) {
    for (gcs::Vertex vertex : mesh.vertices()) {
      positions.push_back(toPoint((*result.vertexPositions)[vertex]));
    }
    return positions;
  }

  if (result.extrinsicGeometry) {
    auto& geometry = const_cast<gcs::VertexPositionGeometry&>(*result.extrinsicGeometry);
    for (gcs::Vertex vertex : mesh.vertices()) {
      positions.push_back(toPoint(geometry.vertexPositions[vertex]));
    }
    return positions;
  }

  throw std::runtime_error("result has no drawable 3D vertex positions; use extrinsic or subdivision retrieval");
}

std::vector<double> vertexScalarData(const gcs::VertexData<double>& field) {
  auto& mesh = *field.getMesh();
  std::vector<double> values;
  values.reserve(mesh.nVertices());
  for (gcs::Vertex vertex : mesh.vertices()) {
    values.push_back(field[vertex]);
  }
  return values;
}

std::vector<double> faceScalarData(const gcs::FaceData<double>& field) {
  auto& mesh = *field.getMesh();
  std::vector<double> values;
  values.reserve(mesh.nFaces());
  for (gcs::Face face : mesh.faces()) {
    values.push_back(field[face]);
  }
  return values;
}

std::vector<Point3> faceVectorData(const gcs::FaceData<geometrycentral::Vector3>& field) {
  auto& mesh = *field.getMesh();
  std::vector<Point3> values;
  values.reserve(mesh.nFaces());
  for (gcs::Face face : mesh.faces()) {
    values.push_back(toPoint(field[face]));
  }
  return values;
}

void addGeneratorTraces(const geodesic_draping::DrapeResult& result) {
  if (!result.extrinsicGenerators) {
    return;
  }

  for (size_t i = 0; i < result.extrinsicGenerators->size(); ++i) {
    const auto& trace = (*result.extrinsicGenerators)[i];
    if (trace.points.size() < 2) {
      continue;
    }
    std::vector<Point3> points;
    points.reserve(trace.points.size());
    for (const geometrycentral::Vector3& point : trace.points) {
      points.push_back(toPoint(point));
    }
    polyscope::registerCurveNetworkLine("generator " + std::to_string(i), points);
  }

  if (result.extrinsicSeed) {
    const std::vector<Point3> seed{toPoint(*result.extrinsicSeed)};
    polyscope::registerPointCloud("seed", seed);
  }
}

void addResultQuantities(polyscope::SurfaceMesh* mesh,
                         const geodesic_draping::DrapeResult& result) {
  if (result.faceShear) {
    mesh->addFaceScalarQuantity("face shear", faceScalarData(*result.faceShear))->setEnabled(true);
  }
  if (result.vertexShear) {
    mesh->addVertexScalarQuantity("vertex shear", vertexScalarData(*result.vertexShear));
  }
  if (result.distances) {
    mesh->addVertexScalarQuantity("distance 0", vertexScalarData((*result.distances)[0]));
    mesh->addVertexScalarQuantity("distance 1", vertexScalarData((*result.distances)[1]));
  }
  if (result.directionFields) {
    mesh->addFaceVectorQuantity("direction field 0", faceVectorData((*result.directionFields)[0]));
    mesh->addFaceVectorQuantity("direction field 1", faceVectorData((*result.directionFields)[1]));
  }
}

} // namespace

int main(int argc, char** argv) {
  try {
    std::string fixtureName = "demo_part";
    geodesic_draping::DrapeSolveOptions solveOptions;
    geodesic_draping::RetrievalOptions retrievalOptions;
    geodesic_draping::IntrinsicConstructionOptions intrinsicOptions;
    geodesic_draping::RefinementOptions refinementOptions;

    solveOptions.mode = geodesic_draping::DrapeSolveMode::Complete;
    retrievalOptions.domain = geodesic_draping::RetrievalDomain::Subdivision;

    for (int i = 1; i < argc; ++i) {
      const std::string arg = argv[i];
      if (arg == "--help" || arg == "-h") {
        printUsage(argv[0]);
        return 0;
      }
      if (arg == "--mode") {
        if (++i >= argc) throw std::runtime_error("--mode requires a value");
        solveOptions.mode = parseMode(argv[i]);
      } else if (arg == "--domain") {
        if (++i >= argc) throw std::runtime_error("--domain requires a value");
        retrievalOptions.domain = parseDomain(argv[i]);
      } else if (arg == "--backend") {
        if (++i >= argc) throw std::runtime_error("--backend requires a value");
        intrinsicOptions.backend = parseBackend(argv[i]);
      } else if (arg == "--refinement") {
        if (++i >= argc) throw std::runtime_error("--refinement requires a value");
        refinementOptions.mode = parseRefinement(argv[i]);
      } else if (arg == "--sample-vertex-shear") {
        retrievalOptions.sampleVertexShear = true;
      } else if (arg.rfind("--", 0) == 0) {
        throw std::runtime_error("unknown option: " + arg);
      } else {
        fixtureName = arg;
      }
    }

    const std::filesystem::path fixtureDir =
        std::filesystem::path(GEODESIC_DRAPING_TEST_DATA_DIR) / fixtureName;
    const geodesic_draping::SurfaceMeshData meshData =
        geodesic_draping::fixture_io::loadMesh(fixtureDir);
    const geodesic_draping::Vec2 seedXY =
        geodesic_draping::fixture_io::loadSeedXY(fixtureDir);
    const double fabricAngle =
        geodesic_draping::fixture_io::loadAngleDegrees(fixtureDir);

    geodesic_draping::DrapeResult result = geodesic_draping::solveDrape(
        meshData,
        seedXY,
        fabricAngle,
        {},
        solveOptions,
        retrievalOptions,
        intrinsicOptions,
        refinementOptions);

    polyscope::init();
    polyscope::options::programName = "Geodesic Draping Result Debug";
    polyscope::SurfaceMesh* psMesh = polyscope::registerSurfaceMesh(
        "drape result",
        vertexPositions(result),
        faceIndices(result));
    addResultQuantities(psMesh, result);
    addGeneratorTraces(result);

    std::cout << "fixture: " << fixtureName << "\n"
              << "vertices: " << result.mesh->nVertices()
              << "  faces: " << result.mesh->nFaces() << "\n";
    polyscope::show();
  } catch (const std::exception& e) {
    std::cerr << "debug_drape_result failed: " << e.what() << "\n\n";
    printUsage(argv[0]);
    return 1;
  }
  return 0;
}

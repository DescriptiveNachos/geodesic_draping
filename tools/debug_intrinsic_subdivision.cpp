#include "fixture_io.h"
#include "geodesic_draping/geometrycentral_adapter.h"

#include "geometrycentral/surface/common_subdivision.h"
#include "geometrycentral/surface/signpost_intrinsic_triangulation.h"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <limits>
#include <optional>
#include <stdexcept>
#include <string>

namespace {

enum class RefinementMode {
  None,
  Flip,
  Refine,
};

struct Options {
  std::string fixtureName = "demo_part";
  RefinementMode refinement = RefinementMode::None;
  std::optional<double> angleThreshold;
  std::optional<double> circumradiusThreshold;
  std::optional<size_t> maxInsertions;
  std::optional<std::string> writeObjPrefix;
};

struct CountStats {
  size_t count = 0;
  size_t nonFinite = 0;
  double min = std::numeric_limits<double>::infinity();
  double max = -std::numeric_limits<double>::infinity();
  double mean = 0.0;
};

void printUsage(const char* argv0) {
  std::cout << "Usage:\n"
            << "  " << argv0
            << " [fixture-name] [--refinement none|flip|refine]\n"
            << "    [--angle-threshold value] [--circumradius-threshold value]\n"
            << "    [--max-insertions value] [--write-obj prefix]\n\n"
            << "Default fixture: demo_part\n"
            << "Fixtures are loaded from " << GEODESIC_DRAPING_TEST_DATA_DIR << "\n";
}

RefinementMode parseRefinementMode(const std::string& value) {
  if (value == "none") {
    return RefinementMode::None;
  }
  if (value == "flip") {
    return RefinementMode::Flip;
  }
  if (value == "refine") {
    return RefinementMode::Refine;
  }
  throw std::runtime_error("--refinement must be one of: none, flip, refine");
}

std::string refinementModeName(RefinementMode mode) {
  switch (mode) {
  case RefinementMode::None:
    return "none";
  case RefinementMode::Flip:
    return "flip";
  case RefinementMode::Refine:
    return "refine";
  }
  return "unknown";
}

Options parseOptions(int argc, char** argv) {
  Options options;
  for (int i = 1; i < argc; ++i) {
    const std::string arg = argv[i];
    if (arg == "--help" || arg == "-h") {
      printUsage(argv[0]);
      std::exit(0);
    }
    if (arg == "--refinement") {
      if (i + 1 >= argc) {
        throw std::runtime_error("--refinement requires a value");
      }
      options.refinement = parseRefinementMode(argv[++i]);
      continue;
    }
    if (arg == "--angle-threshold") {
      if (i + 1 >= argc) {
        throw std::runtime_error("--angle-threshold requires a value");
      }
      options.angleThreshold = std::stod(argv[++i]);
      continue;
    }
    if (arg == "--circumradius-threshold") {
      if (i + 1 >= argc) {
        throw std::runtime_error("--circumradius-threshold requires a value");
      }
      options.circumradiusThreshold = std::stod(argv[++i]);
      continue;
    }
    if (arg == "--max-insertions") {
      if (i + 1 >= argc) {
        throw std::runtime_error("--max-insertions requires a value");
      }
      options.maxInsertions = static_cast<size_t>(std::stoull(argv[++i]));
      continue;
    }
    if (arg == "--write-obj") {
      if (i + 1 >= argc) {
        throw std::runtime_error("--write-obj requires a prefix");
      }
      options.writeObjPrefix = argv[++i];
      continue;
    }
    options.fixtureName = arg;
  }
  if ((options.angleThreshold || options.circumradiusThreshold || options.maxInsertions) &&
      options.refinement != RefinementMode::Refine) {
    throw std::runtime_error("refinement thresholds are only valid with --refinement refine");
  }
  return options;
}

CountStats scalarStats(const std::vector<double>& values) {
  CountStats stats;
  stats.count = values.size();
  double sum = 0.0;
  size_t finiteCount = 0;
  for (double value : values) {
    if (!std::isfinite(value)) {
      ++stats.nonFinite;
      continue;
    }
    stats.min = std::min(stats.min, value);
    stats.max = std::max(stats.max, value);
    sum += value;
    ++finiteCount;
  }
  stats.mean = finiteCount > 0 ? sum / static_cast<double>(finiteCount)
                               : std::numeric_limits<double>::quiet_NaN();
  if (finiteCount == 0) {
    stats.min = std::numeric_limits<double>::quiet_NaN();
    stats.max = std::numeric_limits<double>::quiet_NaN();
  }
  return stats;
}

void printStats(const std::string& label, const CountStats& stats) {
  std::cout << "  " << std::left << std::setw(18) << label
            << " n " << stats.count
            << "  min " << stats.min
            << "  max " << stats.max
            << "  mean " << stats.mean
            << "  nonfinite " << stats.nonFinite << "\n";
}

template <typename MeshT>
size_t countBoundaryEdges(MeshT& mesh) {
  size_t count = 0;
  for (auto e : mesh.edges()) {
    if (e.isBoundary()) {
      ++count;
    }
  }
  return count;
}

void printSubdivisionPointTypeCounts(geometrycentral::surface::CommonSubdivision& subdivision) {
  size_t vertexVertex = 0;
  size_t edgeTransverse = 0;
  size_t edgeParallel = 0;
  size_t faceVertex = 0;
  size_t edgeVertex = 0;
  for (const auto& point : subdivision.subdivisionPoints) {
    switch (point.intersectionType) {
    case geometrycentral::surface::CSIntersectionType::VERTEX_VERTEX:
      ++vertexVertex;
      break;
    case geometrycentral::surface::CSIntersectionType::EDGE_TRANSVERSE:
      ++edgeTransverse;
      break;
    case geometrycentral::surface::CSIntersectionType::EDGE_PARALLEL:
      ++edgeParallel;
      break;
    case geometrycentral::surface::CSIntersectionType::FACE_VERTEX:
      ++faceVertex;
      break;
    case geometrycentral::surface::CSIntersectionType::EDGE_VERTEX:
      ++edgeVertex;
      break;
    }
  }
  std::cout << "common subdivision raw points " << subdivision.subdivisionPoints.size() << "\n"
            << "  point_types"
            << " vertex_vertex " << vertexVertex
            << " edge_transverse " << edgeTransverse
            << " edge_parallel " << edgeParallel
            << " face_vertex " << faceVertex
            << " edge_vertex " << edgeVertex << "\n";
}

template <typename EdgeDataT>
void printPointsAlongStats(const std::string& label, const EdgeDataT& pointsAlong) {
  size_t empty = 0;
  size_t shortPath = 0;
  size_t longest = 0;
  size_t total = 0;
  auto* mesh = pointsAlong.getMesh();
  for (auto e : mesh->edges()) {
    const auto& points = pointsAlong[e];
    if (points.empty()) {
      ++empty;
    }
    if (points.size() < 2) {
      ++shortPath;
    }
    longest = std::max(longest, points.size());
    total += points.size();
  }
  const double mean = mesh->nEdges() > 0
                          ? static_cast<double>(total) / static_cast<double>(mesh->nEdges())
                          : 0.0;
  std::cout << "  " << label
            << " empty " << empty
            << " short " << shortPath
            << " longest " << longest
            << " mean " << mean << "\n";
}

void ensureParentDirectory(const std::string& prefix) {
  const std::filesystem::path path(prefix);
  if (path.has_parent_path()) {
    std::filesystem::create_directories(path.parent_path());
  }
}

} // namespace

int main(int argc, char** argv) {
  try {
    const Options options = parseOptions(argc, argv);
    const std::filesystem::path fixtureDir =
        std::filesystem::path(GEODESIC_DRAPING_TEST_DATA_DIR) / options.fixtureName;
    const geodesic_draping::SurfaceMeshData meshData =
        geodesic_draping::fixture_io::loadMesh(fixtureDir);
    auto surface = geodesic_draping::makeGeometryCentralSurface(meshData);

    std::cout << std::setprecision(17);
    std::cout << "fixture " << options.fixtureName << "\n"
              << "input mesh V/E/F "
              << surface.mesh->nVertices() << "/"
              << surface.mesh->nEdges() << "/"
              << surface.mesh->nFaces()
              << "  boundary_edges " << countBoundaryEdges(*surface.mesh) << "\n"
              << "refinement " << refinementModeName(options.refinement);
    if (options.angleThreshold) {
      std::cout << "  angle_threshold " << *options.angleThreshold;
    }
    if (options.circumradiusThreshold) {
      std::cout << "  circumradius_threshold " << *options.circumradiusThreshold;
    }
    if (options.maxInsertions) {
      std::cout << "  max_insertions " << *options.maxInsertions;
    }
    std::cout << "\n";

    geometrycentral::surface::SignpostIntrinsicTriangulation triangulation(
        *surface.mesh, *surface.geometry);

    if (options.refinement == RefinementMode::Flip) {
      triangulation.flipToDelaunay();
    } else if (options.refinement == RefinementMode::Refine) {
      const double angleThreshold = options.angleThreshold.value_or(25.0);
      const double circumradiusThreshold =
          options.circumradiusThreshold.value_or(std::numeric_limits<double>::infinity());
      const size_t maxInsertions =
          options.maxInsertions.value_or(geometrycentral::INVALID_IND);
      triangulation.delaunayRefine(angleThreshold, circumradiusThreshold, maxInsertions);
    }

    std::cout << "intrinsic mesh V/E/F "
              << triangulation.intrinsicMesh->nVertices() << "/"
              << triangulation.intrinsicMesh->nEdges() << "/"
              << triangulation.intrinsicMesh->nFaces()
              << "  boundary_edges " << countBoundaryEdges(*triangulation.intrinsicMesh) << "\n";

    geometrycentral::surface::CommonSubdivision& subdivision =
        triangulation.getCommonSubdivision();
    printSubdivisionPointTypeCounts(subdivision);
    const auto [expectedV, expectedE, expectedF] = subdivision.elementCounts();
    std::cout << "  expected constructed V/E/F "
              << expectedV << "/" << expectedE << "/" << expectedF << "\n";
    printPointsAlongStats("pointsAlongA", subdivision.pointsAlongA);
    printPointsAlongStats("pointsAlongB", subdivision.pointsAlongB);

    bool simpleMeshBuilt = false;
    try {
      std::unique_ptr<geometrycentral::surface::SimplePolygonMesh> simpleMesh =
          subdivision.buildSimpleMesh();
      simpleMeshBuilt = true;
      std::cout << "buildSimpleMesh ok"
                << "  V/F " << simpleMesh->nVertices()
                << "/" << simpleMesh->nFaces() << "\n";
      if (options.writeObjPrefix) {
        ensureParentDirectory(*options.writeObjPrefix);
        simpleMesh->writeMesh(*options.writeObjPrefix + "_simple.obj", "obj");
      }
    } catch (const std::exception& e) {
      std::cout << "buildSimpleMesh failed  error: " << e.what() << "\n";
    }

    bool constructed = false;
    try {
      subdivision.constructMesh();
      constructed = true;
      std::cout << "constructMesh ok"
                << "  V/E/F " << subdivision.mesh->nVertices()
                << "/" << subdivision.mesh->nEdges()
                << "/" << subdivision.mesh->nFaces()
                << "  boundary_edges " << countBoundaryEdges(*subdivision.mesh) << "\n";
    } catch (const std::exception& e) {
      std::cout << "constructMesh failed  error: " << e.what() << "\n";
    }

    if (constructed) {
      surface.geometry->requireVertexPositions();
      surface.geometry->requireVertexNormals();

      geometrycentral::surface::VertexData<geometrycentral::Vector3> transferredPositions =
          subdivision.interpolateAcrossA(surface.geometry->vertexPositions);
      geometrycentral::surface::VertexData<geometrycentral::Vector3> transferredNormals =
          subdivision.interpolateAcrossA(surface.geometry->vertexNormals);

      std::vector<double> positionNorms;
      std::vector<double> normalNorms;
      positionNorms.reserve(subdivision.mesh->nVertices());
      normalNorms.reserve(subdivision.mesh->nVertices());
      for (auto v : subdivision.mesh->vertices()) {
        positionNorms.push_back(transferredPositions[v].norm());
        normalNorms.push_back(transferredNormals[v].norm());
      }
      std::cout << "transferred input vertex data to subdivision\n";
      printStats("position_norm", scalarStats(positionNorms));
      printStats("normal_norm", scalarStats(normalNorms));

      geometrycentral::surface::FaceData<double> intrinsicColors =
          geometrycentral::surface::niceColors(*triangulation.intrinsicMesh);
      geometrycentral::surface::FaceData<double> subdivisionColors =
          subdivision.copyFromB(intrinsicColors);
      std::vector<double> colors;
      colors.reserve(subdivision.mesh->nFaces());
      for (auto f : subdivision.mesh->faces()) {
        colors.push_back(subdivisionColors[f]);
      }
      printStats("intrinsic_color", scalarStats(colors));

      if (options.writeObjPrefix) {
        ensureParentDirectory(*options.writeObjPrefix);
        subdivision.writeToFile(*options.writeObjPrefix, surface.geometry->vertexPositions);
      }
    } else if (simpleMeshBuilt) {
      std::cout << "data transfer skipped because constructMesh failed; simple mesh still built\n";
    }

    return constructed ? 0 : 2;
  } catch (const std::exception& e) {
    std::cerr << "error: " << e.what() << "\n\n";
    printUsage(argv[0]);
    return 1;
  }
}

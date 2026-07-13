#include "geodesic_draping/geodrape.h"

#include "geodrape_internal.h"

#include "geometrycentral/surface/common_subdivision.h"
#include "geometrycentral/surface/transfer_functions.h"

#include <cmath>
#include <stdexcept>

namespace geodesic_draping {
namespace {
namespace gcs = geometrycentral::surface;

geometrycentral::Vector3 zeroVector() {
  return geometrycentral::Vector3{0.0, 0.0, 0.0};
}

geometrycentral::Vector3 normalized(const geometrycentral::Vector3& value) {
  const double magnitude = geometrycentral::norm(value);
  if (magnitude <= 0.0) {
    return zeroVector();
  }
  return value / magnitude;
}

gcs::FaceData<geometrycentral::Vector3> toFaceData(gcs::SurfaceMesh& mesh,
                                                   const FaceHeatDirectionField& values) {
  if (values.size() != mesh.nFaces()) {
    throw std::runtime_error("face vector field size does not match target mesh");
  }
  gcs::FaceData<geometrycentral::Vector3> data(mesh);
  for (gcs::Face face : mesh.faces()) {
    data[face] = values[face.getIndex()];
  }
  return data;
}

gcs::FaceData<double> toFaceData(gcs::SurfaceMesh& mesh, const std::vector<double>& values) {
  if (values.size() != mesh.nFaces()) {
    throw std::runtime_error("face scalar field size does not match target mesh");
  }
  gcs::FaceData<double> data(mesh);
  for (gcs::Face face : mesh.faces()) {
    data[face] = values[face.getIndex()];
  }
  return data;
}

gcs::VertexData<double> toVertexData(gcs::SurfaceMesh& mesh, const std::vector<double>& values) {
  if (values.size() != mesh.nVertices()) {
    throw std::runtime_error("vertex scalar field size does not match target mesh");
  }
  gcs::VertexData<double> data(mesh);
  for (gcs::Vertex vertex : mesh.vertices()) {
    data[vertex] = values[vertex.getIndex()];
  }
  return data;
}

gcs::VertexData<double> intrinsicFaceDataToVertexData(gcs::SurfaceMesh& mesh,
                                                      const gcs::FaceData<double>& faceData) {
  std::vector<double> faceScalars(mesh.nFaces(), 0.0);
  for (gcs::Face face : mesh.faces()) {
    faceScalars[face.getIndex()] = faceData[face];
  }
  return toVertexData(mesh, averageIntrinsicFaceScalarsToVertices(mesh, faceScalars));
}

geometrycentral::Vector3 inputPosition(const gcs::SurfacePoint& inputPoint,
                                       gcs::VertexPositionGeometry& inputGeometry) {
  return inputPoint.interpolate(inputGeometry.vertexPositions);
}

geometrycentral::Vector3 intrinsicPointToInputPosition(gcs::IntrinsicTriangulation& triangulation,
                                                       gcs::VertexPositionGeometry& inputGeometry,
                                                       const gcs::SurfacePoint& intrinsicPoint) {
  return inputPosition(triangulation.equivalentPointOnInput(intrinsicPoint), inputGeometry);
}

ExtrinsicGeneratorTrace toExtrinsicTrace(gcs::IntrinsicTriangulation& triangulation,
                                         gcs::VertexPositionGeometry& inputGeometry,
                                         const IntrinsicGeneratorTrace& trace) {
  ExtrinsicGeneratorTrace converted;
  converted.hitBoundary = trace.hitBoundary;
  converted.length = trace.length;
  converted.points.reserve(trace.points.size());
  for (const gcs::SurfacePoint& point : trace.points) {
    converted.points.push_back(intrinsicPointToInputPosition(triangulation, inputGeometry, point));
  }
  return converted;
}

std::array<ExtrinsicGeneratorTrace, 4> toExtrinsicTraces(
    gcs::IntrinsicTriangulation& triangulation,
    gcs::VertexPositionGeometry& inputGeometry,
    const std::array<IntrinsicGeneratorTrace, 4>& traces) {
  std::array<ExtrinsicGeneratorTrace, 4> converted;
  for (size_t i = 0; i < traces.size(); ++i) {
    converted[i] = toExtrinsicTrace(triangulation, inputGeometry, traces[i]);
  }
  return converted;
}

std::array<geometrycentral::Vector3, 4> directionsFromExtrinsicTraces(
    const geometrycentral::Vector3& seed,
    const std::array<ExtrinsicGeneratorTrace, 4>& traces) {
  std::array<geometrycentral::Vector3, 4> directions;
  directions.fill(zeroVector());
  for (size_t i = 0; i < traces.size(); ++i) {
    for (const geometrycentral::Vector3& point : traces[i].points) {
      const geometrycentral::Vector3 candidate = point - seed;
      if (geometrycentral::norm(candidate) > 0.0) {
        directions[i] = normalized(candidate);
        break;
      }
    }
  }
  return directions;
}

std::array<gcs::VertexData<double>, 2> transferIntrinsicVertexDataToInput(
    gcs::IntrinsicTriangulation& triangulation,
    const std::array<gcs::VertexData<double>, 2>& values) {
  return {
      gcs::transferBtoA(triangulation, values[0], gcs::TransferMethod::Pointwise),
      gcs::transferBtoA(triangulation, values[1], gcs::TransferMethod::Pointwise),
  };
}

gcs::VertexData<double> transferIntrinsicVertexDataToInput(gcs::IntrinsicTriangulation& triangulation,
                                                           const gcs::VertexData<double>& values) {
  return gcs::transferBtoA(triangulation, values, gcs::TransferMethod::Pointwise);
}

gcs::VertexData<geometrycentral::Vector3> subdivisionVertexPositions(
    gcs::CommonSubdivision& subdivision,
    gcs::VertexPositionGeometry& inputGeometry) {
  subdivision.constructMesh();
  gcs::VertexData<geometrycentral::Vector3> positions(*subdivision.mesh);
  for (gcs::Vertex vertex : subdivision.mesh->vertices()) {
    positions[vertex] = subdivision.sourcePoints[vertex]->posA.interpolate(inputGeometry.vertexPositions);
  }
  return positions;
}

} // namespace

IntrinsicDrapeResult GeoDrapeSolver::retrieveIntrinsic(bool sampleVertexShear) const {
  const CoreIntrinsicResult& core = lastResult();
  gcs::SurfaceMesh& mesh = *intrinsicTriangulation_->intrinsicMesh;

  IntrinsicDrapeResult result;
  result.mesh = &mesh;
  result.geometry = intrinsicTriangulation_.get();
  result.mode = core.mode;
  result.seed = core.intrinsicSeed;
  result.directions = core.intrinsicDirections;
  result.generators = core.generators;
  result.directionFields = {
      toFaceData(mesh, core.directions[0]),
      toFaceData(mesh, core.directions[1]),
  };
  if (core.distances) {
    result.distances = std::array<gcs::VertexData<double>, 2>{
        toVertexData(mesh, (*core.distances)[0]),
        toVertexData(mesh, (*core.distances)[1]),
    };
  }
  if (core.faceShear) {
    result.faceShear = toFaceData(mesh, *core.faceShear);
    if (sampleVertexShear) {
      result.vertexShear = intrinsicFaceDataToVertexData(mesh, *result.faceShear);
    }
  }
  return result;
}

ExtrinsicDrapeResult GeoDrapeSolver::retrieveExtrinsic(bool sampleVertexShear) const {
  const CoreIntrinsicResult& core = lastResult();
  gcs::SurfaceMesh& intrinsicMesh = *intrinsicTriangulation_->intrinsicMesh;

  ExtrinsicDrapeResult result;
  result.mesh = inputSurface_.mesh.get();
  result.geometry = inputSurface_.geometry.get();
  result.mode = core.mode;
  result.seed = intrinsicPointToInputPosition(
      *intrinsicTriangulation_,
      *inputSurface_.geometry,
      core.intrinsicSeed);
  result.generators = toExtrinsicTraces(
      *intrinsicTriangulation_,
      *inputSurface_.geometry,
      core.generators);
  result.directions = directionsFromExtrinsicTraces(result.seed, result.generators);

  if (core.distances) {
    const std::array<gcs::VertexData<double>, 2> intrinsicDistances{
        toVertexData(intrinsicMesh, (*core.distances)[0]),
        toVertexData(intrinsicMesh, (*core.distances)[1]),
    };
    result.distances = transferIntrinsicVertexDataToInput(*intrinsicTriangulation_, intrinsicDistances);
  }
  if (sampleVertexShear && core.faceShear) {
    const gcs::FaceData<double> intrinsicFaceShear = toFaceData(intrinsicMesh, *core.faceShear);
    const gcs::VertexData<double> intrinsicVertexShear =
        intrinsicFaceDataToVertexData(intrinsicMesh, intrinsicFaceShear);
    result.vertexShear =
        transferIntrinsicVertexDataToInput(*intrinsicTriangulation_, intrinsicVertexShear);
  }
  return result;
}

SubdivisionDrapeResult GeoDrapeSolver::retrieveSubdivision(bool sampleVertexShear) const {
  const CoreIntrinsicResult& core = lastResult();
  gcs::SurfaceMesh& intrinsicMesh = *intrinsicTriangulation_->intrinsicMesh;
  gcs::CommonSubdivision& subdivision = intrinsicTriangulation_->getCommonSubdivision();
  subdivision.constructMesh();

  SubdivisionDrapeResult result;
  result.mesh = subdivision.mesh.get();
  result.mode = core.mode;
  result.vertexPositions = subdivisionVertexPositions(subdivision, *inputSurface_.geometry);
  result.seed = intrinsicPointToInputPosition(
      *intrinsicTriangulation_,
      *inputSurface_.geometry,
      core.intrinsicSeed);
  result.generators = toExtrinsicTraces(
      *intrinsicTriangulation_,
      *inputSurface_.geometry,
      core.generators);
  result.directions = directionsFromExtrinsicTraces(result.seed, result.generators);

  const gcs::FaceData<geometrycentral::Vector3> intrinsicDirections0 =
      toFaceData(intrinsicMesh, core.directions[0]);
  const gcs::FaceData<geometrycentral::Vector3> intrinsicDirections1 =
      toFaceData(intrinsicMesh, core.directions[1]);
  result.directionFields = {
      subdivision.copyFromB(intrinsicDirections0),
      subdivision.copyFromB(intrinsicDirections1),
  };

  if (core.distances) {
    const std::array<gcs::VertexData<double>, 2> intrinsicDistances{
        toVertexData(intrinsicMesh, (*core.distances)[0]),
        toVertexData(intrinsicMesh, (*core.distances)[1]),
    };
    result.distances = std::array<gcs::VertexData<double>, 2>{
        subdivision.interpolateAcrossB(intrinsicDistances[0]),
        subdivision.interpolateAcrossB(intrinsicDistances[1]),
    };
  }
  if (core.faceShear) {
    const gcs::FaceData<double> intrinsicFaceShear = toFaceData(intrinsicMesh, *core.faceShear);
    result.faceShear = subdivision.copyFromB(intrinsicFaceShear);
    if (sampleVertexShear) {
      result.vertexShear = intrinsicFaceDataToVertexData(*subdivision.mesh, *result.faceShear);
    }
  }
  return result;
}

} // namespace geodesic_draping

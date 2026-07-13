#include "geodesic_draping/geodrape.h"

#include "geodrape_internal.h"

#include "geometrycentral/surface/common_subdivision.h"
#include "geometrycentral/surface/transfer_functions.h"

#include <Eigen/Dense>

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

geometrycentral::Vector3 barycentricDisplacementToSubdivisionVector(
    gcs::CommonSubdivision& subdivision,
    gcs::VertexPositionGeometry& inputGeometry,
    gcs::Face subdivisionFace,
    gcs::Face intrinsicFace,
    const geometrycentral::Vector3& barycentricDisplacement) {
  Eigen::Matrix3d barycentricBasis;
  Eigen::Matrix<double, 3, 3> positionBasis;

  size_t i = 0;
  for (gcs::Vertex vertex : subdivisionFace.adjacentVertices()) {
    if (i >= 3) {
      throw std::runtime_error("subdivision vector transfer expects triangular faces");
    }
    const gcs::SurfacePoint pointOnIntrinsic =
        subdivision.sourcePoints[vertex]->posB.inFace(intrinsicFace);
    const geometrycentral::Vector3 pointOnInput =
        subdivision.sourcePoints[vertex]->posA.interpolate(inputGeometry.vertexPositions);

    barycentricBasis(0, i) = pointOnIntrinsic.faceCoords.x;
    barycentricBasis(1, i) = pointOnIntrinsic.faceCoords.y;
    barycentricBasis(2, i) = pointOnIntrinsic.faceCoords.z;

    positionBasis(0, i) = pointOnInput.x;
    positionBasis(1, i) = pointOnInput.y;
    positionBasis(2, i) = pointOnInput.z;
    ++i;
  }
  if (i != 3) {
    throw std::runtime_error("subdivision vector transfer expects triangular faces");
  }

  const Eigen::Vector3d displacement(
      barycentricDisplacement.x,
      barycentricDisplacement.y,
      barycentricDisplacement.z);
  // The heat field stores a barycentric displacement on the intrinsic source
  // face; map that displacement through the local affine common-subdivision
  // face before exposing it as a drawable 3D vector.
  const Eigen::Vector3d subdivisionWeights =
      barycentricBasis.colPivHouseholderQr().solve(displacement);
  const Eigen::Vector3d vector = positionBasis * subdivisionWeights;
  return normalized(geometrycentral::Vector3{vector.x(), vector.y(), vector.z()});
}

gcs::FaceData<geometrycentral::Vector3> subdivisionDirectionField(
    gcs::CommonSubdivision& subdivision,
    gcs::VertexPositionGeometry& inputGeometry,
    const FaceHeatDirectionField& intrinsicDirections) {
  gcs::FaceData<geometrycentral::Vector3> field(*subdivision.mesh);
  for (gcs::Face face : subdivision.mesh->faces()) {
    const gcs::Face intrinsicFace = subdivision.sourceFaceB[face];
    field[face] = barycentricDisplacementToSubdivisionVector(
        subdivision,
        inputGeometry,
        face,
        intrinsicFace,
        intrinsicDirections[intrinsicFace.getIndex()]);
  }
  return field;
}

} // namespace

DrapeResult GeoDrapeSolver::retrieveIntrinsic(bool sampleVertexShear) const {
  const CoreIntrinsicResult& core = lastCoreResult();
  gcs::SurfaceMesh& mesh = *intrinsicTriangulation_->intrinsicMesh;

  DrapeResult result;
  result.domain = RetrievalDomain::Intrinsic;
  result.mesh = &mesh;
  result.intrinsicGeometry = intrinsicTriangulation_.get();
  result.mode = core.mode;
  result.intrinsicSeed = core.intrinsicSeed;
  result.intrinsicDirections = core.intrinsicDirections;
  result.intrinsicGenerators = core.generators;
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

DrapeResult GeoDrapeSolver::retrieveExtrinsic(bool sampleVertexShear) const {
  const CoreIntrinsicResult& core = lastCoreResult();
  gcs::SurfaceMesh& intrinsicMesh = *intrinsicTriangulation_->intrinsicMesh;

  DrapeResult result;
  result.domain = RetrievalDomain::Extrinsic;
  result.mesh = inputSurface_.mesh.get();
  result.extrinsicGeometry = inputSurface_.geometry.get();
  result.mode = core.mode;
  result.extrinsicSeed = intrinsicPointToInputPosition(
      *intrinsicTriangulation_,
      *inputSurface_.geometry,
      core.intrinsicSeed);
  result.extrinsicGenerators = toExtrinsicTraces(
      *intrinsicTriangulation_,
      *inputSurface_.geometry,
      core.generators);
  result.extrinsicDirections = directionsFromExtrinsicTraces(*result.extrinsicSeed, *result.extrinsicGenerators);

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

DrapeResult GeoDrapeSolver::retrieveSubdivision(bool sampleVertexShear) const {
  const CoreIntrinsicResult& core = lastCoreResult();
  gcs::SurfaceMesh& intrinsicMesh = *intrinsicTriangulation_->intrinsicMesh;
  gcs::CommonSubdivision& subdivision = intrinsicTriangulation_->getCommonSubdivision();
  subdivision.constructMesh();

  DrapeResult result;
  result.domain = RetrievalDomain::Subdivision;
  result.mesh = subdivision.mesh.get();
  result.mode = core.mode;
  result.vertexPositions = subdivisionVertexPositions(subdivision, *inputSurface_.geometry);
  result.extrinsicSeed = intrinsicPointToInputPosition(
      *intrinsicTriangulation_,
      *inputSurface_.geometry,
      core.intrinsicSeed);
  result.extrinsicGenerators = toExtrinsicTraces(
      *intrinsicTriangulation_,
      *inputSurface_.geometry,
      core.generators);
  result.extrinsicDirections = directionsFromExtrinsicTraces(*result.extrinsicSeed, *result.extrinsicGenerators);

  result.directionFields = {
      subdivisionDirectionField(subdivision, *inputSurface_.geometry, core.directions[0]),
      subdivisionDirectionField(subdivision, *inputSurface_.geometry, core.directions[1]),
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

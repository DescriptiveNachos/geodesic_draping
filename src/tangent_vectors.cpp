#include "geodrape_internal.h"

#include "geometrycentral/surface/trace_geodesic.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace geodesic_draping {
namespace {

geometrycentral::Vector3 toGcVector3(const Vec3& value) {
  return geometrycentral::Vector3{value.x(), value.y(), value.z()};
}

gcs::BarycentricVector embeddedActiveFaceDirection(gcs::VertexPositionGeometry& inputGeometry,
                                                  gcs::IntrinsicTriangulation& triangulation,
                                                  gcs::Face activeFace,
                                                  double angle) {
  std::array<Vec3, 3> positions;
  size_t i = 0;
  for (gcs::Vertex vertex : activeFace.adjacentVertices()) {
    const gcs::SurfacePoint inputPoint = triangulation.equivalentPointOnInput(gcs::SurfacePoint(vertex));
    positions[i] = interpolateSurfacePoint(inputPoint, inputGeometry);
    ++i;
  }

  const Vec3 ambient(std::cos(angle * kPi / 180.0),
                     std::sin(angle * kPi / 180.0),
                     0.0);
  const Vec3 e1 = positions[1] - positions[0];
  const Vec3 e2 = positions[2] - positions[0];
  Eigen::Matrix2d gram;
  gram << e1.dot(e1), e1.dot(e2),
          e2.dot(e1), e2.dot(e2);
  Eigen::Vector2d rhs;
  rhs << ambient.dot(e1), ambient.dot(e2);
  const Eigen::Vector2d uv = gram.ldlt().solve(rhs);
  return normalizeVector(
      gcs::BarycentricVector(activeFace, geometrycentral::Vector3{-uv.x() - uv.y(), uv.x(), uv.y()}),
      triangulation);
}

} // namespace

TraceSettings resolveTraceSettings(const TraceSettings& defaults,
                                   const AdvancedTraceOptions& overrides) {
  TraceSettings resolved = defaults;
  if (overrides.traceLength) {
    resolved.traceLength = *overrides.traceLength;
  }
  if (overrides.maxIterations) {
    resolved.maxIterations = *overrides.maxIterations;
  }
  return resolved;
}

double boundingBoxDiagonal(const SurfaceMeshData& meshData) {
  if (meshData.vertices.empty()) {
    return 1.0;
  }

  Vec3 minCorner = meshData.vertices.front();
  Vec3 maxCorner = meshData.vertices.front();
  for (const Vec3& vertex : meshData.vertices) {
    minCorner = minCorner.cwiseMin(vertex);
    maxCorner = maxCorner.cwiseMax(vertex);
  }
  const double diagonal = (maxCorner - minCorner).norm();
  return diagonal > 0.0 ? diagonal : 1.0;
}

TraceSettings makeTraceDefaults(const SurfaceMeshData& meshData,
                                const gcs::SurfaceMesh& activeMesh) {
  TraceSettings settings;
  settings.traceLength = 100.0 * boundingBoxDiagonal(meshData);
  settings.maxIterations = std::max<size_t>(1, 100 * activeMesh.nFaces());
  return settings;
}

gcs::SurfacePoint toFaceSurfacePoint(gcs::SurfaceMesh& mesh, const BarycentricPoint& point) {
  return gcs::SurfacePoint(
      mesh.face(point.faceIndex),
      geometrycentral::Vector3{
          point.barycentric.x(),
          point.barycentric.y(),
          point.barycentric.z(),
      });
}

gcs::BarycentricVector normalizeVector(gcs::BarycentricVector vector,
                                       gcs::IntrinsicGeometryInterface& geometry) {
  const double magnitude = vector.norm(geometry);
  if (magnitude == 0.0) {
    throw std::runtime_error("cannot normalize zero tangent direction");
  }
  return vector / magnitude;
}

gcs::BarycentricVector intrinsicDirectionFromFabricAngle(const SurfaceMeshData& meshData,
                                                        gcs::VertexPositionGeometry& inputGeometry,
                                                        gcs::IntrinsicTriangulation& triangulation,
                                                        const gcs::SurfacePoint& inputSeed,
                                                        const gcs::SurfacePoint& intrinsicSeed,
                                                        double fabricAngle,
                                                        bool inputConnectivityPreserved,
                                                        bool useCommonSubdivisionInputAdapter) {
  const gcs::SurfacePoint inputFaceSeed = inputSeed.inSomeFace();
  inputGeometry.requireFaceTangentBasis();

  const Vec3 ambientDirection(std::cos(fabricAngle * kPi / 180.0),
                              std::sin(fabricAngle * kPi / 180.0),
                              0.0);
  const geometrycentral::Vector3 ambient = toGcVector3(ambientDirection);
  const geometrycentral::Vector3 basisX = inputGeometry.faceTangentBasis[inputFaceSeed.face][0];
  const geometrycentral::Vector3 basisY = inputGeometry.faceTangentBasis[inputFaceSeed.face][1];
  geometrycentral::Vector2 tangentDirection{
      dot(ambient, basisX),
      dot(ambient, basisY),
  };
  inputGeometry.unrequireFaceTangentBasis();

  const double tangentNorm = norm(tangentDirection);
  if (tangentNorm == 0.0) {
    throw std::runtime_error("fabric direction is degenerate in the seed tangent plane");
  }
  tangentDirection /= tangentNorm;

  const gcs::SurfacePoint intrinsicSeedFace = intrinsicSeed.inSomeFace();
  if (inputConnectivityPreserved) {
    return embeddedActiveFaceDirection(inputGeometry, triangulation, intrinsicSeedFace.face, fabricAngle);
  }

  const double baseStep = std::max(1e-9, 1e-6 * boundingBoxDiagonal(meshData));
  gcs::TraceOptions options;
  options.includePath = false;
  options.errorOnProblem = false;
  options.maxIters = 16;

  for (size_t attempt = 0; attempt < 8; ++attempt) {
    const double step = baseStep * std::pow(0.25, static_cast<double>(attempt));
    const gcs::TraceGeodesicResult inputTrace =
        gcs::traceGeodesic(inputGeometry, inputSeed, step * tangentDirection, options);
    const gcs::SurfacePoint intrinsicEndpoint =
        inputToIntrinsic(triangulation, inputTrace.endPoint, useCommonSubdivisionInputAdapter).inSomeFace();
    try {
      const gcs::SurfacePoint endpointInSeedFace = intrinsicEndpoint.inFace(intrinsicSeedFace.face);
      return normalizeVector(gcs::BarycentricVector(intrinsicSeedFace, endpointInSeedFace),
                             triangulation);
    } catch (const std::exception&) {
      // Retry with a shorter local displacement if the mapped endpoint crossed a face boundary.
    }
  }

  return embeddedActiveFaceDirection(inputGeometry, triangulation, intrinsicSeedFace.face, fabricAngle);
}

} // namespace geodesic_draping

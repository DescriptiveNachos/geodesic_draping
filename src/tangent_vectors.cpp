#include "geodrape_internal.h"

#include "geodesic_draping/signed_heat.h"

#include "geometrycentral/surface/trace_geodesic.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace geodesic_draping {
namespace {

geometrycentral::Vector3 toGcVector3(const Vec3& value) {
  return geometrycentral::Vector3{value.x(), value.y(), value.z()};
}

gcs::BarycentricVector embeddedActiveFaceDirection(ReferenceGeometry& reference,
                                                  ActiveIntrinsicDomain& activeDomain,
                                                  gcs::Face activeFace,
                                                  double angle) {
  std::array<Vec3, 3> positions;
  size_t i = 0;
  for (gcs::Vertex vertex : activeFace.adjacentVertices()) {
    const gcs::SurfacePoint inputPoint = activeDomain.intrinsicToInput(gcs::SurfacePoint(vertex));
    positions[i] = interpolateSurfacePoint(inputPoint, *reference.surface().geometry);
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
      activeDomain.geometry());
}

GeneratorTrace traceActiveGenerator(ReferenceGeometry& reference,
                                    ActiveIntrinsicDomain& activeDomain,
                                    const gcs::SurfacePoint& start,
                                    const gcs::BarycentricVector& direction,
                                    const TraceSettings& settings) {
  gcs::SurfacePoint startFace = start.inSomeFace();
  gcs::BarycentricVector directionInFace = direction.inFace(startFace.face);
  directionInFace = normalizeVector(directionInFace, activeDomain.geometry()) * settings.traceLength;

  gcs::TraceOptions options;
  options.includePath = true;
  options.errorOnProblem = false;
  options.barrierEdges = nullptr;
  options.maxIters = settings.maxIterations;

  const gcs::TraceGeodesicResult traced = gcs::traceGeodesic(
      activeDomain.geometry(),
      startFace.face,
      startFace.faceCoords,
      directionInFace.faceCoords,
      options);
  if (!traced.hasPath) {
    throw std::runtime_error("active traceGeodesic did not return a path");
  }

  GeneratorTrace trace;
  trace.hitBoundary = traced.hitBoundary;
  trace.length = traced.length;
  trace.points.reserve(traced.pathPoints.size());
  trace.surfaceReferences.reserve(traced.pathPoints.size());
  for (const gcs::SurfacePoint& intrinsicPoint : traced.pathPoints) {
    trace.surfaceReferences.push_back(toSurfaceReference(intrinsicPoint));
    const gcs::SurfacePoint inputPoint = activeDomain.intrinsicToInput(intrinsicPoint);
    trace.points.push_back(interpolateSurfacePoint(inputPoint, *reference.surface().geometry));
  }
  return trace;
}

} // namespace

std::array<Vec3, 4> generateCartesianFamilyDirections(double fabricAngle,
                                                       double fiberAngle) {
  const double fabricAngleRadians = fabricAngle * kPi / 180.0;
  const double fiberAngleRadians = fiberAngle * kPi / 180.0;

  const Vec3 direction0(std::cos(fabricAngleRadians), std::sin(fabricAngleRadians), 0.0);
  const Vec3 direction1(std::cos(fabricAngleRadians + fiberAngleRadians),
                        std::sin(fabricAngleRadians + fiberAngleRadians),
                        0.0);
  return {direction0, -direction0, direction1, -direction1};
}

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

gcs::BarycentricVector toBarycentricVector(gcs::SurfaceMesh& mesh,
                                           const TangentVectorRef& ref) {
  if (ref.type != SurfaceReferenceType::Face || ref.coords.size() != 3) {
    throw std::runtime_error("solveFromIntrinsic currently requires a face tangent vector with three barycentric coordinates");
  }
  return gcs::BarycentricVector(
      mesh.face(ref.elementIndex),
      geometrycentral::Vector3{ref.coords[0], ref.coords[1], ref.coords[2]});
}

TangentVectorRef toTangentVectorRef(const gcs::BarycentricVector& vector) {
  const gcs::BarycentricVector faceVector = vector.inSomeFace();
  TangentVectorRef ref;
  ref.type = SurfaceReferenceType::Face;
  ref.elementIndex = faceVector.face.getIndex();
  ref.coords = {faceVector.faceCoords.x, faceVector.faceCoords.y, faceVector.faceCoords.z};
  return ref;
}

std::array<Vec3, 4> cartesianDirectionsFromIntrinsic(
    ReferenceGeometry& reference,
    ActiveIntrinsicDomain& activeDomain,
    const gcs::SurfacePoint& seed,
    const std::array<gcs::BarycentricVector, 4>& directions) {
  std::array<Vec3, 4> out{};
  const gcs::SurfacePoint seedFace = seed.inSomeFace();
  std::array<Vec3, 3> positions;
  size_t vertexIndex = 0;
  for (gcs::Vertex vertex : seedFace.face.adjacentVertices()) {
    const gcs::SurfacePoint inputPoint = activeDomain.intrinsicToInput(gcs::SurfacePoint(vertex));
    positions[vertexIndex] = interpolateSurfacePoint(inputPoint, *reference.surface().geometry);
    ++vertexIndex;
  }

  for (size_t i = 0; i < directions.size(); ++i) {
    const gcs::BarycentricVector inFace = directions[i].inFace(seedFace.face);
    Vec3 delta = inFace.faceCoords.x * positions[0] +
                 inFace.faceCoords.y * positions[1] +
                 inFace.faceCoords.z * positions[2];
    const double norm = delta.norm();
    out[i] = Vec3::Zero();
    if (norm > 0.0) {
      out[i] = delta / norm;
    }
  }
  return out;
}

gcs::BarycentricVector intrinsicDirectionFromFabricAngle(ReferenceGeometry& reference,
                                                        ActiveIntrinsicDomain& activeDomain,
                                                        const gcs::SurfacePoint& inputSeed,
                                                        const gcs::SurfacePoint& intrinsicSeed,
                                                        double fabricAngle,
                                                        bool inputConnectivityPreserved) {
  auto& inputGeometry = *reference.surface().geometry;
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
    return embeddedActiveFaceDirection(reference, activeDomain, intrinsicSeedFace.face, fabricAngle);
  }

  const double baseStep = std::max(1e-9, 1e-6 * boundingBoxDiagonal(reference.meshData()));
  gcs::TraceOptions options;
  options.includePath = false;
  options.errorOnProblem = false;
  options.maxIters = 16;

  for (size_t attempt = 0; attempt < 8; ++attempt) {
    const double step = baseStep * std::pow(0.25, static_cast<double>(attempt));
    const gcs::TraceGeodesicResult inputTrace =
        gcs::traceGeodesic(inputGeometry, inputSeed, step * tangentDirection, options);
    const gcs::SurfacePoint intrinsicEndpoint =
        activeDomain.inputToIntrinsic(inputTrace.endPoint).inSomeFace();
    try {
      const gcs::SurfacePoint endpointInSeedFace = intrinsicEndpoint.inFace(intrinsicSeedFace.face);
      return normalizeVector(gcs::BarycentricVector(intrinsicSeedFace, endpointInSeedFace),
                             activeDomain.geometry());
    } catch (const std::exception&) {
      // Retry with a shorter local displacement if the mapped endpoint crossed a face boundary.
    }
  }

  return embeddedActiveFaceDirection(reference, activeDomain, intrinsicSeedFace.face, fabricAngle);
}

std::array<GeneratorTrace, 4> traceActiveGenerators(
    ReferenceGeometry& reference,
    ActiveIntrinsicDomain& activeDomain,
    const gcs::SurfacePoint& start,
    const std::array<gcs::BarycentricVector, 4>& directions,
    const TraceSettings& settings) {
  return {
      traceActiveGenerator(reference, activeDomain, start, directions[0], settings),
      traceActiveGenerator(reference, activeDomain, start, directions[1], settings),
      traceActiveGenerator(reference, activeDomain, start, directions[2], settings),
      traceActiveGenerator(reference, activeDomain, start, directions[3], settings),
  };
}

} // namespace geodesic_draping

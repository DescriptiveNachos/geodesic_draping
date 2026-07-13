#include "geodrape_internal.h"

#include "geometrycentral/surface/trace_geodesic.h"

#include <stdexcept>

namespace geodesic_draping {
namespace {

geometrycentral::Vector3 toGcVector3(const Vec3& value) {
  return geometrycentral::Vector3{value.x(), value.y(), value.z()};
}

Vec3 fromGcVector3(const geometrycentral::Vector3& value) {
  return Vec3(value.x, value.y, value.z);
}

} // namespace

Vec3 interpolateSurfacePoint(const geometrycentral::surface::SurfacePoint& point,
                             geometrycentral::surface::VertexPositionGeometry& geometry) {
  geometry.requireVertexPositions();
  return fromGcVector3(point.interpolate(geometry.vertexPositions));
}

SurfaceReference toSurfaceReference(const geometrycentral::surface::SurfacePoint& point) {
  SurfaceReference ref;
  switch (point.type) {
  case geometrycentral::surface::SurfacePointType::Vertex:
    ref.type = SurfaceReferenceType::Vertex;
    ref.elementIndex = point.vertex.getIndex();
    break;
  case geometrycentral::surface::SurfacePointType::Edge:
    ref.type = SurfaceReferenceType::Edge;
    ref.elementIndex = point.edge.getIndex();
    ref.params = {point.tEdge};
    break;
  case geometrycentral::surface::SurfacePointType::Face:
    ref.type = SurfaceReferenceType::Face;
    ref.elementIndex = point.face.getIndex();
    ref.params = {point.faceCoords.x, point.faceCoords.y, point.faceCoords.z};
    break;
  }
  return ref;
}

GeneratorTrace traceGeneratorFromFace(GeometryCentralSurface& surface,
                                      const BarycentricPoint& start,
                                      const Vec3& direction,
                                      const TraceSettings& settings) {
  if (!surface.mesh || !surface.geometry) {
    throw std::runtime_error("traceGeneratorFromFace requires a valid geometry-central surface");
  }
  if (start.faceIndex >= surface.mesh->nFaces()) {
    throw std::runtime_error("traceGeneratorFromFace start face index is out of range");
  }

  auto& mesh = *surface.mesh;
  auto& geometry = *surface.geometry;

  geometry.requireFaceTangentBasis();

  const geometrycentral::surface::Face face = mesh.face(start.faceIndex);
  const geometrycentral::Vector3 bary{
      start.barycentric.x(),
      start.barycentric.y(),
      start.barycentric.z(),
  };
  const geometrycentral::Vector3 ambientDirection = settings.traceLength * toGcVector3(direction);
  const geometrycentral::Vector3 basisX = geometry.faceTangentBasis[face][0];
  const geometrycentral::Vector3 basisY = geometry.faceTangentBasis[face][1];
  const geometrycentral::Vector2 tangentDirection{
      dot(ambientDirection, basisX),
      dot(ambientDirection, basisY),
  };

  geometrycentral::surface::TraceOptions options;
  options.includePath = true;
  options.errorOnProblem = false;
  options.barrierEdges = nullptr;
  options.maxIters = settings.maxIterations;

  const geometrycentral::surface::SurfacePoint startPoint(face, bary);
  const geometrycentral::surface::TraceGeodesicResult result =
      geometrycentral::surface::traceGeodesic(geometry, startPoint, tangentDirection, options);

  if (!result.hasPath) {
    throw std::runtime_error("geometry-central traceGeodesic did not return a path");
  }

  GeneratorTrace trace;
  trace.hitBoundary = result.hitBoundary;
  trace.length = result.length;
  trace.points.reserve(result.pathPoints.size());
  trace.surfaceReferences.reserve(result.pathPoints.size());
  for (const geometrycentral::surface::SurfacePoint& point : result.pathPoints) {
    trace.points.push_back(interpolateSurfacePoint(point, geometry));
    trace.surfaceReferences.push_back(toSurfaceReference(point));
  }
  return trace;
}

std::array<GeneratorTrace, 4> traceGenerators(GeometryCentralSurface& surface,
                                              const BarycentricPoint& start,
                                              const std::array<Vec3, 4>& directions,
                                              const TraceSettings& settings) {
  return {
      traceGeneratorFromFace(surface, start, directions[0], settings),
      traceGeneratorFromFace(surface, start, directions[1], settings),
      traceGeneratorFromFace(surface, start, directions[2], settings),
      traceGeneratorFromFace(surface, start, directions[3], settings),
  };
}

IntrinsicGeneratorTrace traceIntrinsicGenerator(gcs::IntrinsicTriangulation& triangulation,
                                                const gcs::SurfacePoint& start,
                                                const gcs::BarycentricVector& direction,
                                                const TraceSettings& settings) {
  gcs::SurfacePoint startFace = start.inSomeFace();
  gcs::BarycentricVector directionInFace = direction.inFace(startFace.face);
  directionInFace = normalizeVector(directionInFace, triangulation) * settings.traceLength;

  gcs::TraceOptions options;
  options.includePath = true;
  options.errorOnProblem = false;
  options.barrierEdges = nullptr;
  options.maxIters = settings.maxIterations;

  const gcs::TraceGeodesicResult traced = gcs::traceGeodesic(
      triangulation,
      startFace.face,
      startFace.faceCoords,
      directionInFace.faceCoords,
      options);
  if (!traced.hasPath) {
    throw std::runtime_error("intrinsic traceGeodesic did not return a path");
  }

  IntrinsicGeneratorTrace trace;
  trace.hitBoundary = traced.hitBoundary;
  trace.length = traced.length;
  trace.points = traced.pathPoints;
  return trace;
}

std::array<IntrinsicGeneratorTrace, 4> traceIntrinsicGenerators(
    gcs::IntrinsicTriangulation& triangulation,
    const gcs::SurfacePoint& start,
    const std::array<gcs::BarycentricVector, 4>& directions,
    const TraceSettings& settings) {
  return {
      traceIntrinsicGenerator(triangulation, start, directions[0], settings),
      traceIntrinsicGenerator(triangulation, start, directions[1], settings),
      traceIntrinsicGenerator(triangulation, start, directions[2], settings),
      traceIntrinsicGenerator(triangulation, start, directions[3], settings),
  };
}

std::array<gcs::Curve, 2> pairOppositeIntrinsicGeneratorTraces(
    const std::array<IntrinsicGeneratorTrace, 4>& traces) {
  std::array<gcs::Curve, 2> curves;
  for (size_t pairIndex = 0; pairIndex < curves.size(); ++pairIndex) {
    const IntrinsicGeneratorTrace& negativeTrace = traces[2 * pairIndex + 1];
    const IntrinsicGeneratorTrace& positiveTrace = traces[2 * pairIndex];
    gcs::Curve& curve = curves[pairIndex];
    curve.isSigned = true;
    curve.nodes.reserve(negativeTrace.points.size() + positiveTrace.points.size());
    curve.nodes.insert(curve.nodes.end(), negativeTrace.points.rbegin(), negativeTrace.points.rend());
    curve.nodes.insert(curve.nodes.end(), positiveTrace.points.begin(), positiveTrace.points.end());
  }
  return curves;
}

SourceCurves pairOppositeGeneratorTraces(const std::array<GeneratorTrace, 4>& traces) {
  SourceCurves sourceCurves;

  for (size_t pairIndex = 0; pairIndex < 2; ++pairIndex) {
    const GeneratorTrace& negativeTrace = traces[2 * pairIndex + 1];
    const GeneratorTrace& positiveTrace = traces[2 * pairIndex];
    auto& curve = sourceCurves.curves[pairIndex];
    curve.reserve(negativeTrace.surfaceReferences.size() + positiveTrace.surfaceReferences.size());
    curve.insert(curve.end(), negativeTrace.surfaceReferences.rbegin(), negativeTrace.surfaceReferences.rend());
    curve.insert(curve.end(), positiveTrace.surfaceReferences.begin(), positiveTrace.surfaceReferences.end());
  }

  return sourceCurves;
}

} // namespace geodesic_draping

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

} // namespace geodesic_draping

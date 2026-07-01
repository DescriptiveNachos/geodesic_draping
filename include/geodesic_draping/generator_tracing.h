#pragma once

#include "geodesic_draping/geometrycentral_adapter.h"
#include "geodesic_draping/seed_projection.h"

#include "geometrycentral/surface/surface_point.h"

#include <array>
#include <optional>
#include <vector>

namespace geodesic_draping {

enum class SurfaceReferenceType {
  Vertex = 0,
  Edge = 1,
  Face = 2,
};

struct SurfaceReference {
  SurfaceReferenceType type = SurfaceReferenceType::Vertex;
  size_t elementIndex = 0;
  std::vector<double> params;
};

struct GeneratorTrace {
  std::vector<Vec3> points;
  std::vector<SurfaceReference> surfaceReferences;
  bool hitBoundary = false;
  double length = 0.0;
};

struct TraceSettings {
  double traceLength = 10000.0;
  size_t maxIterations = geometrycentral::INVALID_IND;
};

GeneratorTrace traceGeneratorFromFace(GeometryCentralSurface& surface,
                                      const BarycentricPoint& start,
                                      const Vec3& direction,
                                      const TraceSettings& settings = {});

std::array<GeneratorTrace, 4> traceGenerators(GeometryCentralSurface& surface,
                                              const BarycentricPoint& start,
                                              const std::array<Vec3, 4>& directions,
                                              const TraceSettings& settings = {});

Vec3 interpolateSurfacePoint(const geometrycentral::surface::SurfacePoint& point,
                             geometrycentral::surface::VertexPositionGeometry& geometry);

SurfaceReference toSurfaceReference(const geometrycentral::surface::SurfacePoint& point);

} // namespace geodesic_draping

#pragma once

#include "geodesic_draping/generator_tracing.h"

#include "geometrycentral/surface/signed_heat_method.h"

#include <string>
#include <vector>

namespace geodesic_draping {

struct SignedHeatSolveOptions {
  bool preserveSourceNormals = false;
  geometrycentral::LevelSetConstraint levelSetConstraint = geometrycentral::LevelSetConstraint::None;
  double softLevelSetWeight = -1.0;
  double diffusionTimeCoefficient = 1.0;
};

geometrycentral::surface::SurfacePoint toGeometryCentralSurfacePoint(
    geometrycentral::surface::SurfaceMesh& mesh,
    const SurfaceReference& ref);

geometrycentral::surface::Curve toGeometryCentralCurve(
    geometrycentral::surface::SurfaceMesh& mesh,
    const std::vector<SurfaceReference>& refs,
    bool isSigned = true);

} // namespace geodesic_draping

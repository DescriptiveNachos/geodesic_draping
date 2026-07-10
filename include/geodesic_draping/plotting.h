#pragma once

#include "geodesic_draping/generator_tracing.h"
#include "geodesic_draping/geodrape.h"
#include "geodesic_draping/seed_projection.h"

#include <array>
#include <string>

namespace geodesic_draping {

struct ProjectionPlotOptions {
  std::string name = "geodesic draping seed";
  double directionLength = 25.0;
  bool clearExisting = false;
  bool show = true;
};

void plotSeedProjectionStep(const SurfaceMeshData& mesh,
                            const SeedProjection& projection,
                            const std::array<Vec3, 4>& directions,
                            const ProjectionPlotOptions& options = {});

void plotGeneratorTraces(const std::array<GeneratorTrace, 4>& traces,
                         const ProjectionPlotOptions& options = {});

void plotDrapeComparisonResult(const SurfaceMeshData& mesh,
                               const DrapeResult& completeResult,
                               const DrapeResult& fastResult,
                               const DrapeResult& overlayResult,
                               const ProjectionPlotOptions& options = {});

} // namespace geodesic_draping

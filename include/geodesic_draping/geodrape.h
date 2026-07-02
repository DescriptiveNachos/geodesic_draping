#pragma once

#include "geodesic_draping/field_processing.h"
#include "geodesic_draping/generator_tracing.h"
#include "geodesic_draping/signed_heat.h"
#include "geodesic_draping/signed_vector_heat.h"

#include <array>
#include <vector>

namespace geodesic_draping {

struct CompleteDrapeResult {
  SeedProjection seed;
  std::array<Vec3, 4> directions;
  std::array<GeneratorTrace, 4> generators;
  SourceCurves sourceCurves;
  std::array<std::vector<double>, 2> distances;
  std::array<std::vector<Vec3>, 2> gradients;
  std::vector<double> shearAnglesDegrees;
};

struct FastDrapeResult {
  SeedProjection seed;
  std::array<Vec3, 4> directions;
  std::array<GeneratorTrace, 4> generators;
  SourceCurves sourceCurves;
  std::array<CustomSignedHeatResult, 2> customHeatSolves;
  std::array<FaceHeatDirectionField, 2> faceDirections;
  std::array<std::vector<double>, 2> distances;
  std::vector<double> faceShearAnglesDegrees;
  std::vector<double> vertexShearAnglesDegrees;
};

struct FastDrapeOptions {
  bool returnDistances = false;
};

class GeoDrapeSolver {
public:
  explicit GeoDrapeSolver(SurfaceMeshData meshData,
                          const SignedHeatSolveOptions& heatOptions = {});

  CompleteDrapeResult solveComplete(const Vec2& seedXY, double angleDegrees);
  FastDrapeResult solveFast(const Vec2& seedXY,
                            double angleDegrees,
                            const FastDrapeOptions& fastOptions = {});

private:
  SurfaceMeshData meshData_;
  GeometryCentralSurface surface_;
  SignedHeatSolveOptions heatOptions_;
  SignedHeatDistanceSolver distanceSolver_;
  CustomSignedHeatSolver customHeatSolver_;
};

CompleteDrapeResult solveCompleteDrape(const SurfaceMeshData& mesh,
                                       const Vec2& seedXY,
                                       double angleDegrees,
                                       const SignedHeatSolveOptions& heatOptions = {});

FastDrapeResult solveFastDrape(const SurfaceMeshData& mesh,
                               const Vec2& seedXY,
                               double angleDegrees,
                               const SignedHeatSolveOptions& heatOptions = {},
                               const FastDrapeOptions& fastOptions = {});

} // namespace geodesic_draping

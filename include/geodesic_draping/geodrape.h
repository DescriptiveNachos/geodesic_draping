#pragma once

#include "geodesic_draping/field_processing.h"
#include "geodesic_draping/generator_tracing.h"
#include "geodesic_draping/custom_signed_heat.h"

#include <array>
#include <optional>
#include <vector>

namespace geodesic_draping {

enum class DrapeSolveMode {
  Fast,
  Hybrid,
  Complete,
};

struct DrapeSolveOptions {
  DrapeSolveMode mode = DrapeSolveMode::Fast;
  bool sampleSecondaryShear = false;
};

struct DrapeResult {
  DrapeSolveMode mode = DrapeSolveMode::Fast;
  SeedProjection seed;
  std::array<Vec3, 4> directions;
  std::array<GeneratorTrace, 4> generators;
  SourceCurves sourceCurves;
  std::array<CustomSignedHeatResult, 2> customHeatSolves;
  std::array<FaceHeatDirectionField, 2> faceDirections;
  std::optional<std::array<std::vector<double>, 2>> distances;
  std::optional<std::array<std::vector<Vec3>, 2>> gradients;
  std::optional<std::vector<double>> faceShearAnglesDegrees;
  std::optional<std::vector<double>> vertexShearAnglesDegrees;
};

class GeoDrapeSolver {
public:
  explicit GeoDrapeSolver(SurfaceMeshData meshData,
                          const SignedHeatSolveOptions& heatOptions = {});

  DrapeResult solve(const Vec2& seedXY,
                    double angleDegrees,
                    const DrapeSolveOptions& solveOptions = {});

private:
  SurfaceMeshData meshData_;
  GeometryCentralSurface surface_;
  SignedHeatSolveOptions heatOptions_;
  CustomSignedHeatSolver customHeatSolver_;
};

DrapeResult solveDrape(const SurfaceMeshData& mesh,
                       const Vec2& seedXY,
                       double angleDegrees,
                       const SignedHeatSolveOptions& heatOptions = {},
                       const DrapeSolveOptions& solveOptions = {});

} // namespace geodesic_draping

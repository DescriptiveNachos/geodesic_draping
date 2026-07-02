#include "geodesic_draping/geodrape.h"

#include <stdexcept>
#include <utility>

namespace geodesic_draping {

GeoDrapeSolver::GeoDrapeSolver(SurfaceMeshData meshData,
                               const SignedHeatSolveOptions& heatOptions)
    : meshData_(std::move(meshData)),
      surface_(makeGeometryCentralSurface(meshData_)),
      heatOptions_(heatOptions),
      customHeatSolver_(surface_, heatOptions_.diffusionTimeCoefficient) {}

DrapeResult GeoDrapeSolver::solve(const Vec2& seedXY,
                                  double angleDegrees,
                                  const DrapeSolveOptions& solveOptions) {
  DrapeResult result;

  const std::optional<SeedProjection> seed = projectPointXYToMesh(meshData_, seedXY);
  if (!seed) {
    throw std::runtime_error("GeoDrapeSolver failed to project seed point to mesh");
  }
  result.seed = *seed;
  result.directions = generateOrthogonalDirections(angleDegrees);

  result.generators = traceGenerators(surface_, result.seed.surfacePoint, result.directions);
  result.sourceCurves = pairOppositeGeneratorTraces(result.generators);
  const bool computeDistances = solveOptions.mode != DrapeSolveMode::Fast;
  result.customHeatSolves = customHeatSolver_.solve(result.sourceCurves, heatOptions_, computeDistances);
  for (size_t i = 0; i < result.customHeatSolves.size(); ++i) {
    result.faceDirections[i] = result.customHeatSolves[i].normalizedFaceDirections;
    if (computeDistances) {
      result.distances[i] = result.customHeatSolves[i].distance;
    }
  }
  result.faceShearAnglesDegrees =
      computeFaceShearAnglesDegrees(surface_, result.faceDirections[0], result.faceDirections[1]);
  result.vertexShearAnglesDegrees =
      averageFaceScalarsToVertices(meshData_, result.faceShearAnglesDegrees, FaceScalarAveraging::FaceArea);
  if (solveOptions.mode == DrapeSolveMode::Complete) {
    result.gradients[0] = computeVertexScalarGradients(meshData_, result.distances[0]);
    result.gradients[1] = computeVertexScalarGradients(meshData_, result.distances[1]);
    result.vertexShearAnglesDegrees = computeShearAnglesDegrees(result.gradients[0], result.gradients[1]);
  }

  return result;
}

CompleteDrapeResult GeoDrapeSolver::solveComplete(const Vec2& seedXY, double angleDegrees) {
  const auto unified = solve(seedXY, angleDegrees, {DrapeSolveMode::Complete});

  CompleteDrapeResult result;
  result.seed = unified.seed;
  result.directions = unified.directions;
  result.generators = unified.generators;
  result.sourceCurves = unified.sourceCurves;
  result.distances = unified.distances;
  result.gradients = unified.gradients;
  result.shearAnglesDegrees = unified.vertexShearAnglesDegrees;

  return result;
}

FastDrapeResult GeoDrapeSolver::solveFast(const Vec2& seedXY,
                                          double angleDegrees,
                                          const DrapeSolveOptions& solveOptions) {
  const auto unified = solve(seedXY, angleDegrees, solveOptions);

  FastDrapeResult result;
  result.seed = unified.seed;
  result.directions = unified.directions;
  result.generators = unified.generators;
  result.sourceCurves = unified.sourceCurves;
  result.customHeatSolves = unified.customHeatSolves;
  result.faceDirections = unified.faceDirections;
  result.distances = unified.distances;
  result.faceShearAnglesDegrees = unified.faceShearAnglesDegrees;
  result.vertexShearAnglesDegrees = unified.vertexShearAnglesDegrees;
  return result;
}

CompleteDrapeResult solveCompleteDrape(const SurfaceMeshData& mesh,
                                       const Vec2& seedXY,
                                       double angleDegrees,
                                       const SignedHeatSolveOptions& heatOptions) {
  GeoDrapeSolver solver(mesh, heatOptions);
  return solver.solveComplete(seedXY, angleDegrees);
}

FastDrapeResult solveFastDrape(const SurfaceMeshData& mesh,
                               const Vec2& seedXY,
                               double angleDegrees,
                               const SignedHeatSolveOptions& heatOptions,
                               const DrapeSolveOptions& solveOptions) {
  GeoDrapeSolver solver(mesh, heatOptions);
  return solver.solveFast(seedXY, angleDegrees, solveOptions);
}

} // namespace geodesic_draping

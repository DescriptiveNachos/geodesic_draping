#include "geodesic_draping/geodrape.h"

#include <stdexcept>
#include <utility>

namespace geodesic_draping {

GeoDrapeSolver::GeoDrapeSolver(SurfaceMeshData meshData,
                               const SignedHeatSolveOptions& heatOptions)
    : meshData_(std::move(meshData)),
      surface_(makeGeometryCentralSurface(meshData_)),
      heatOptions_(heatOptions),
      distanceSolver_(surface_, heatOptions_.diffusionTimeCoefficient),
      customHeatSolver_(surface_, heatOptions_.diffusionTimeCoefficient) {}

CompleteDrapeResult GeoDrapeSolver::solveComplete(const Vec2& seedXY, double angleDegrees) {
  CompleteDrapeResult result;

  const std::optional<SeedProjection> seed = projectPointXYToMesh(meshData_, seedXY);
  if (!seed) {
    throw std::runtime_error("solveCompleteDrape failed to project seed point to mesh");
  }
  result.seed = *seed;
  result.directions = generateOrthogonalDirections(angleDegrees);

  result.generators = traceGenerators(surface_, result.seed.surfacePoint, result.directions);
  result.sourceCurves = pairOppositeGeneratorTraces(result.generators);
  result.distances = distanceSolver_.computeDistances(result.sourceCurves, heatOptions_);

  result.gradients[0] = computeVertexScalarGradients(meshData_, result.distances[0]);
  result.gradients[1] = computeVertexScalarGradients(meshData_, result.distances[1]);
  result.shearAnglesDegrees = computeShearAnglesDegrees(result.gradients[0], result.gradients[1]);

  return result;
}

FastDrapeResult GeoDrapeSolver::solveFast(const Vec2& seedXY,
                                          double angleDegrees,
                                          const FastDrapeOptions& fastOptions) {
  FastDrapeResult result;

  const std::optional<SeedProjection> seed = projectPointXYToMesh(meshData_, seedXY);
  if (!seed) {
    throw std::runtime_error("solveFastDrape failed to project seed point to mesh");
  }
  result.seed = *seed;
  result.directions = generateOrthogonalDirections(angleDegrees);

  result.generators = traceGenerators(surface_, result.seed.surfacePoint, result.directions);
  result.sourceCurves = pairOppositeGeneratorTraces(result.generators);
  for (size_t i = 0; i < result.customHeatSolves.size(); ++i) {
    result.customHeatSolves[i] =
        customHeatSolver_.solve(result.sourceCurves.curves[i], heatOptions_, fastOptions.returnDistances);
    result.faceDirections[i] = result.customHeatSolves[i].normalizedFaceDirections;
    if (fastOptions.returnDistances) {
      result.distances[i] = result.customHeatSolves[i].distance;
    }
  }
  result.faceShearAnglesDegrees =
      computeFaceShearAnglesDegrees(surface_, result.faceDirections[0], result.faceDirections[1]);
  result.vertexShearAnglesDegrees =
      averageFaceScalarsToVertices(meshData_, result.faceShearAnglesDegrees, FaceScalarAveraging::FaceArea);

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
                               const FastDrapeOptions& fastOptions) {
  GeoDrapeSolver solver(mesh, heatOptions);
  return solver.solveFast(seedXY, angleDegrees, fastOptions);
}

} // namespace geodesic_draping

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

FastDrapeResult GeoDrapeSolver::solveFast(const Vec2& seedXY, double angleDegrees) {
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
    result.customHeatSolves[i].diffusion =
        customHeatSolver_.solveDiffusedEdgeHeatField(result.sourceCurves.curves[i], heatOptions_);
    result.customHeatSolves[i].normalizedFaceDirections =
        sampleAndNormalizeFaceDirections(surface_, result.customHeatSolves[i].diffusion.diffusedEdgeHeatField);
    result.customHeatSolves[i].vertexDirections =
        averageFaceDirectionsToVerticesReference(surface_, result.customHeatSolves[i].normalizedFaceDirections);
  }
  result.gradients[0] = result.customHeatSolves[0].vertexDirections;
  result.gradients[1] = result.customHeatSolves[1].vertexDirections;
  result.shearAnglesDegrees = computeShearAnglesDegrees(result.gradients[0], result.gradients[1]);

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
                               const SignedHeatSolveOptions& heatOptions) {
  GeoDrapeSolver solver(mesh, heatOptions);
  return solver.solveFast(seedXY, angleDegrees);
}

} // namespace geodesic_draping

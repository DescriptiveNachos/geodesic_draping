#include "geodesic_draping/geodrape.h"

#include <stdexcept>

namespace geodesic_draping {

CompleteDrapeResult solveCompleteDrape(const SurfaceMeshData& mesh,
                                       const Vec2& seedXY,
                                       double angleDegrees,
                                       const SignedHeatSolveOptions& heatOptions) {
  CompleteDrapeResult result;

  const std::optional<SeedProjection> seed = projectPointXYToMesh(mesh, seedXY);
  if (!seed) {
    throw std::runtime_error("solveCompleteDrape failed to project seed point to mesh");
  }
  result.seed = *seed;
  result.directions = generateOrthogonalDirections(angleDegrees);

  GeometryCentralSurface surface = makeGeometryCentralSurface(mesh);
  result.generators = traceGenerators(surface, result.seed.surfacePoint, result.directions);
  result.sourceCurves = pairOppositeGeneratorTraces(result.generators);
  result.distances = computeSignedHeatDistances(surface, result.sourceCurves, heatOptions);

  result.gradients[0] = computeVertexScalarGradients(mesh, result.distances[0]);
  result.gradients[1] = computeVertexScalarGradients(mesh, result.distances[1]);
  result.shearAnglesDegrees = computeShearAnglesDegrees(result.gradients[0], result.gradients[1]);

  return result;
}

FastDrapeResult solveFastDrape(const SurfaceMeshData& mesh,
                               const Vec2& seedXY,
                               double angleDegrees,
                               const SignedHeatSolveOptions& heatOptions) {
  FastDrapeResult result;

  const std::optional<SeedProjection> seed = projectPointXYToMesh(mesh, seedXY);
  if (!seed) {
    throw std::runtime_error("solveFastDrape failed to project seed point to mesh");
  }
  result.seed = *seed;
  result.directions = generateOrthogonalDirections(angleDegrees);

  GeometryCentralSurface surface = makeGeometryCentralSurface(mesh);
  result.generators = traceGenerators(surface, result.seed.surfacePoint, result.directions);
  result.sourceCurves = pairOppositeGeneratorTraces(result.generators);
  result.heatVectorSolves = computeSignedVectorHeats(surface, result.sourceCurves, heatOptions);
  result.gradients[0] = result.heatVectorSolves[0].vertexVectorHeat;
  result.gradients[1] = result.heatVectorSolves[1].vertexVectorHeat;
  result.shearAnglesDegrees = computeShearAnglesDegrees(result.gradients[0], result.gradients[1]);

  return result;
}

} // namespace geodesic_draping

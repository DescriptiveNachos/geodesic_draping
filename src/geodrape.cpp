#include "geodesic_draping/geodrape.h"

#include "geodrape_internal.h"

#include <stdexcept>

namespace geodesic_draping {

GeoDrapeSolver::GeoDrapeSolver(SurfaceMeshData meshData,
                               const SignedHeatSolveOptions& heatOptions)
    : GeoDrapeSolver(std::move(meshData), heatOptions, {}, {}) {}

GeoDrapeSolver::GeoDrapeSolver(SurfaceMeshData meshData,
                               const SignedHeatSolveOptions& heatOptions,
                               const IntrinsicConstructionOptions& intrinsicOptions,
                               const RefinementOptions& refinementOptions)
    : meshData_(std::move(meshData)),
      inputSurface_(makeGeometryCentralSurface(meshData_)),
      heatOptions_(heatOptions),
      inputConnectivityPreserved_(refinementOptions.mode == RefinementMode::None) {
  useCommonSubdivisionInputAdapter_ =
      intrinsicOptions.backend == IntrinsicTriangulationBackend::IntegerCoordinates;
  intrinsicTriangulation_ = makeIntrinsicTriangulation(
      intrinsicOptions.backend,
      *inputSurface_.mesh,
      *inputSurface_.geometry);
  applyRefinement(*intrinsicTriangulation_, refinementOptions);
  traceDefaults_ = makeTraceDefaults(meshData_, *intrinsicTriangulation_->intrinsicMesh);
  customHeatSolver_ = std::make_unique<CustomSignedHeatSolver>(
      *intrinsicTriangulation_->intrinsicMesh,
      *intrinsicTriangulation_,
      heatOptions_.diffusionTimeCoefficient);
}

const CoreIntrinsicResult& GeoDrapeSolver::solve(const Vec2& seedXY,
                                                 double fabricAngle,
                                                 const DrapeSolveOptions& solveOptions) {
  const TraceSettings trace = resolveTraceSettings(traceDefaults_, solveOptions.advanced.trace);
  const IntrinsicSolveInput input = adaptExtrinsicInput(
      seedXY,
      fabricAngle,
      solveOptions.fiberAngle,
      solveOptions.mode,
      trace);
  lastIntrinsicResult_ = solveCore(input);
  return *lastIntrinsicResult_;
}

const CoreIntrinsicResult& GeoDrapeSolver::solveFromIntrinsic(
    const gcs::SurfacePoint& seed,
    const gcs::BarycentricVector& fabricDirection,
    double fiberAngle,
    const DrapeSolveOptions& solveOptions) {
  const TraceSettings trace = resolveTraceSettings(traceDefaults_, solveOptions.advanced.trace);
  const gcs::SurfacePoint intrinsicSeed = seed.inSomeFace();
  const gcs::BarycentricVector direction0 = normalizeVector(
      fabricDirection.inFace(intrinsicSeed.face),
      *intrinsicTriangulation_);
  gcs::BarycentricVector direction1 =
      normalizeVector(direction0.rotate(*intrinsicTriangulation_, fiberAngle * kPi / 180.0),
                      *intrinsicTriangulation_);

  IntrinsicSolveInput input;
  input.seed = intrinsicSeed;
  input.directions = {direction0, -direction0, direction1, -direction1};
  input.mode = solveOptions.mode;
  input.trace = trace;

  lastIntrinsicResult_ = solveCore(input);
  return *lastIntrinsicResult_;
}

const CoreIntrinsicResult& GeoDrapeSolver::lastResult() const {
  if (!lastIntrinsicResult_) {
    throw std::runtime_error("GeoDrapeSolver::lastResult() requires a previous solve");
  }
  return *lastIntrinsicResult_;
}

IntrinsicSolveInput GeoDrapeSolver::adaptExtrinsicInput(const Vec2& seedXY,
                                                        double fabricAngle,
                                                        double fiberAngle,
                                                        DrapeSolveMode mode,
                                                        const TraceSettings& trace) {

  const std::optional<SeedProjection> seed = projectPointXYToMesh(meshData_, seedXY);
  if (!seed) {
    throw std::runtime_error("GeoDrapeSolver failed to project seed point to mesh");
  }
  const gcs::SurfacePoint inputSeed = toFaceSurfacePoint(*inputSurface_.mesh, seed->surfacePoint);
  const gcs::SurfacePoint intrinsicSeed =
      inputToIntrinsic(*intrinsicTriangulation_, inputSeed, useCommonSubdivisionInputAdapter_).inSomeFace();
  const gcs::BarycentricVector direction0 = intrinsicDirectionFromFabricAngle(
      meshData_,
      *inputSurface_.geometry,
      *intrinsicTriangulation_,
      inputSeed,
      intrinsicSeed,
      fabricAngle,
      inputConnectivityPreserved_,
      useCommonSubdivisionInputAdapter_);
  const gcs::BarycentricVector direction1 = intrinsicDirectionFromFabricAngle(
      meshData_,
      *inputSurface_.geometry,
      *intrinsicTriangulation_,
      inputSeed,
      intrinsicSeed,
      fabricAngle + fiberAngle,
      inputConnectivityPreserved_,
      useCommonSubdivisionInputAdapter_);

  IntrinsicSolveInput input;
  input.seed = intrinsicSeed;
  input.directions = {direction0, -direction0, direction1, -direction1};
  input.mode = mode;
  input.trace = trace;
  return input;
}

CoreIntrinsicResult GeoDrapeSolver::solveCore(const IntrinsicSolveInput& input) {
  CoreIntrinsicResult result;
  result.mode = input.mode;
  result.intrinsicSeed = input.seed;
  result.intrinsicDirections = input.directions;

  result.generators = traceIntrinsicGenerators(
      *intrinsicTriangulation_,
      input.seed,
      input.directions,
      input.trace);
  const std::array<gcs::Curve, 2> sourceCurves = pairOppositeIntrinsicGeneratorTraces(result.generators);

  const bool computeDistances = input.mode != DrapeSolveMode::Fast;
  const std::array<CustomSignedHeatResult, 2> customHeatSolves =
      customHeatSolver_->solve(sourceCurves, heatOptions_, computeDistances);

  if (computeDistances) {
    result.distances = std::array<std::vector<double>, 2>{};
  }
  for (size_t i = 0; i < customHeatSolves.size(); ++i) {
    result.directions[i] = customHeatSolves[i].normalizedFaceDirections;
    if (computeDistances) {
      (*result.distances)[i] = customHeatSolves[i].distance;
    }
  }

  if (input.mode == DrapeSolveMode::Complete) {
    result.directions[0] = computeIntrinsicFaceScalarGradients(
        *intrinsicTriangulation_->intrinsicMesh,
        *intrinsicTriangulation_,
        (*result.distances)[0]);
    result.directions[1] = computeIntrinsicFaceScalarGradients(
        *intrinsicTriangulation_->intrinsicMesh,
        *intrinsicTriangulation_,
        (*result.distances)[1]);
  }
  result.faceShear = computeFaceShearAnglesDegrees(
      *intrinsicTriangulation_->intrinsicMesh,
      *intrinsicTriangulation_,
      result.directions[0],
      result.directions[1]);
  return result;
}

} // namespace geodesic_draping

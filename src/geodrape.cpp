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
    : reference_(std::move(meshData)),
      activeDomain_(reference_, intrinsicOptions, refinementOptions),
      traceDefaults_(makeTraceDefaults(reference_.meshData(), activeDomain_.mesh())),
      heatOptions_(heatOptions),
      customHeatSolver_(std::make_unique<CustomSignedHeatSolver>(
      activeDomain_.mesh(),
      activeDomain_.geometry(),
      heatOptions_.diffusionTimeCoefficient)),
      inputConnectivityPreserved_(refinementOptions.mode == RefinementMode::None) {}

DrapeResult GeoDrapeSolver::solve(const Vec2& seedXY,
                                  double fabricAngle,
                                  const DrapeSolveOptions& solveOptions) {
  return solve(seedXY, fabricAngle, 90.0, solveOptions);
}

DrapeResult GeoDrapeSolver::solve(const Vec2& seedXY,
                                  double fabricAngle,
                                  double fiberAngle,
                                  const DrapeSolveOptions& solveOptions) {
  const TraceSettings trace = resolveTraceSettings(traceDefaults_, solveOptions.advanced.trace);
  const IntrinsicSolveInput input = adaptExtrinsicInput(
      seedXY,
      fabricAngle,
      fiberAngle,
      solveOptions.mode,
      trace);
  lastIntrinsicResult_ = solveCore(input);
  return retrieveFromCore(
      *lastIntrinsicResult_,
      solveOptions.retrieval,
      solveOptions.sampleVertexShear);
}

DrapeResult GeoDrapeSolver::solveFromIntrinsic(const SurfaceReference& seed,
                                               const TangentVectorRef& fabricDirection,
                                               double fiberAngle,
                                               const DrapeSolveOptions& solveOptions) {
  const TraceSettings trace = resolveTraceSettings(traceDefaults_, solveOptions.advanced.trace);
  gcs::SurfacePoint intrinsicSeed =
      toGeometryCentralSurfacePoint(activeDomain_.mesh(), seed).inSomeFace();
  gcs::BarycentricVector direction0 = normalizeVector(
      toBarycentricVector(activeDomain_.mesh(), fabricDirection).inFace(intrinsicSeed.face),
      activeDomain_.geometry());
  gcs::BarycentricVector direction1 =
      normalizeVector(direction0.rotate(activeDomain_.geometry(), fiberAngle * kPi / 180.0),
                      activeDomain_.geometry());

  IntrinsicSolveInput input;
  input.seed = intrinsicSeed;
  input.directions = {direction0, -direction0, direction1, -direction1};
  input.cartesianDirections = cartesianDirectionsFromIntrinsic(
      reference_,
      activeDomain_,
      intrinsicSeed,
      input.directions);
  input.mode = solveOptions.mode;
  input.trace = trace;

  lastIntrinsicResult_ = solveCore(input);
  return retrieveFromCore(*lastIntrinsicResult_, solveOptions.retrieval, solveOptions.sampleVertexShear);
}

DrapeResult GeoDrapeSolver::retrieve(ResultDomain retrieval, bool sampleVertexShear) {
  if (!lastIntrinsicResult_) {
    throw std::runtime_error("GeoDrapeSolver::retrieve() requires a previous solve");
  }
  return retrieveFromCore(*lastIntrinsicResult_, retrieval, sampleVertexShear);
}

IntrinsicSolveInput GeoDrapeSolver::adaptExtrinsicInput(const Vec2& seedXY,
                                                        double fabricAngle,
                                                        double fiberAngle,
                                                        DrapeSolveMode mode,
                                                        const TraceSettings& trace) {
  const std::optional<SeedProjection> seed = projectPointXYToMesh(reference_.meshData(), seedXY);
  if (!seed) {
    throw std::runtime_error("GeoDrapeSolver failed to project seed point to mesh");
  }
  const gcs::SurfacePoint inputSeed = toFaceSurfacePoint(*reference_.surface().mesh, seed->surfacePoint);
  const gcs::SurfacePoint intrinsicSeed = activeDomain_.inputToIntrinsic(inputSeed).inSomeFace();
  const gcs::BarycentricVector direction0 = intrinsicDirectionFromFabricAngle(
      reference_,
      activeDomain_,
      inputSeed,
      intrinsicSeed,
      fabricAngle,
      inputConnectivityPreserved_);
  const gcs::BarycentricVector direction1 = intrinsicDirectionFromFabricAngle(
      reference_,
      activeDomain_,
      inputSeed,
      intrinsicSeed,
      fabricAngle + fiberAngle,
      inputConnectivityPreserved_);

  IntrinsicSolveInput input;
  input.seed = intrinsicSeed;
  input.directions = {direction0, -direction0, direction1, -direction1};
  input.cartesianDirections = generateCartesianFamilyDirections(fabricAngle, fiberAngle);
  input.mode = mode;
  input.trace = trace;
  return input;
}

CoreIntrinsicResult GeoDrapeSolver::solveCore(const IntrinsicSolveInput& input) {
  CoreIntrinsicResult result;
  result.mode = input.mode;
  result.intrinsicSeed = toSurfaceReference(input.seed);
  result.intrinsicDirections = {
      toTangentVectorRef(input.directions[0]),
      toTangentVectorRef(input.directions[1]),
      toTangentVectorRef(input.directions[2]),
      toTangentVectorRef(input.directions[3]),
  };
  result.cartesianDirections = input.cartesianDirections;

  result.generators = traceActiveGenerators(
      reference_,
      activeDomain_,
      input.seed,
      input.directions,
      input.trace);
  const SourceCurves sourceCurves = pairOppositeGeneratorTraces(result.generators);

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
        activeDomain_.mesh(),
        activeDomain_.geometry(),
        (*result.distances)[0]);
    result.directions[1] = computeIntrinsicFaceScalarGradients(
        activeDomain_.mesh(),
        activeDomain_.geometry(),
        (*result.distances)[1]);
  }

  result.faceShear = computeFaceShearAnglesDegrees(
      activeDomain_.mesh(),
      activeDomain_.geometry(),
      result.directions[0],
      result.directions[1]);
  return result;
}

DrapeResult solveDrape(const SurfaceMeshData& mesh,
                       const Vec2& seedXY,
                       double fabricAngle,
                       const SignedHeatSolveOptions& heatOptions,
                       const DrapeSolveOptions& solveOptions) {
  GeoDrapeSolver solver(mesh, heatOptions);
  return solver.solve(seedXY, fabricAngle, solveOptions);
}

} // namespace geodesic_draping

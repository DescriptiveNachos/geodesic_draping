#include "geodesic_draping/geodrape.h"

#include "geodrape_internal.h"

#include <memory>
#include <stdexcept>
#include <utility>

namespace geodesic_draping {

DrapeResult GeoDrapeSolver::solve(const Vec2& seedXY,
                                  double fabricAngle,
                                  const DrapeSolveOptions& solveOptions,
                                  const RetrievalOptions& retrievalOptions) {
  const TraceSettings trace = resolveTraceSettings(traceDefaults_, solveOptions.trace);
  const IntrinsicSolveInput input = adaptExtrinsicInput(
      seedXY,
      fabricAngle,
      solveOptions.fiberAngle,
      solveOptions.mode,
      trace);
  lastIntrinsicResult_ = solveCore(input);
  return retrieve(retrievalOptions);
}

DrapeResult GeoDrapeSolver::solveFromIntrinsic(
    const geometrycentral::surface::SurfacePoint& seed,
    const geometrycentral::surface::BarycentricVector& fabricDirection,
    const DrapeSolveOptions& solveOptions,
    const RetrievalOptions& retrievalOptions) {
  const TraceSettings trace = resolveTraceSettings(traceDefaults_, solveOptions.trace);
  const IntrinsicSolveInput input = adaptIntrinsicInput(seed, fabricDirection, solveOptions, trace);
  lastIntrinsicResult_ = solveCore(input);
  return retrieve(retrievalOptions);
}

DrapeResult GeoDrapeSolver::retrieve(const RetrievalOptions& retrievalOptions) const {
  if (!lastIntrinsicResult_) {
    throw std::runtime_error("GeoDrapeSolver::retrieve() requires a previous solve");
  }

  switch (retrievalOptions.domain) {
  case RetrievalDomain::Intrinsic:
    return retrieveIntrinsic(retrievalOptions.sampleVertexShear);
  case RetrievalDomain::Extrinsic:
    return retrieveExtrinsic(retrievalOptions.sampleVertexShear);
  case RetrievalDomain::Subdivision:
    return retrieveSubdivision(retrievalOptions.sampleVertexShear);
  }

  throw std::runtime_error("GeoDrapeSolver::retrieve() received an unsupported retrieval domain");
}

const CoreIntrinsicResult& GeoDrapeSolver::lastCoreResult() const {
  if (!lastIntrinsicResult_) {
    throw std::runtime_error("GeoDrapeSolver::lastCoreResult() requires a previous solve");
  }
  return *lastIntrinsicResult_;
}

DrapeResult solveDrape(SurfaceMeshData meshData,
                       const Vec2& seedXY,
                       double fabricAngle,
                       const SignedHeatSolveOptions& heatOptions,
                       const DrapeSolveOptions& solveOptions,
                       const RetrievalOptions& retrievalOptions,
                       IntrinsicTriangulationBackend intrinsicBackend,
                       const RefinementOptions& refinementOptions) {
  auto solver = std::make_shared<GeoDrapeSolver>(
      std::move(meshData),
      heatOptions,
      intrinsicBackend,
      refinementOptions);
  DrapeResult result = solver->solve(seedXY, fabricAngle, solveOptions, retrievalOptions);
  result.storageOwner = solver;
  return result;
}

} // namespace geodesic_draping

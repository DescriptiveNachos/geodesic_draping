#include "geodesic_draping/geodrape.h"

#include "geometrycentral/surface/signpost_intrinsic_triangulation.h"

#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>

namespace geodesic_draping {
namespace gcs = geometrycentral::surface;

namespace {

constexpr double pi = 3.141592653589793238462643383279502884;

std::array<Vec3, 4> generateFamilyDirections(double fabricAngleDegrees,
                                             double fiberAngleDegrees) {
  const double fabricAngleRadians = fabricAngleDegrees * pi / 180.0;
  const double fiberAngleRadians = fiberAngleDegrees * pi / 180.0;

  const Vec3 direction0(std::cos(fabricAngleRadians), std::sin(fabricAngleRadians), 0.0);
  const Vec3 direction1(std::cos(fabricAngleRadians + fiberAngleRadians),
                        std::sin(fabricAngleRadians + fiberAngleRadians),
                        0.0);
  return {direction0, -direction0, direction1, -direction1};
}

TraceSettings resolveTraceSettings(const TraceSettings& defaults,
                                   const AdvancedTraceOptions& overrides) {
  TraceSettings resolved = defaults;
  if (overrides.traceLength) {
    resolved.traceLength = *overrides.traceLength;
  }
  if (overrides.maxIterations) {
    resolved.maxIterations = *overrides.maxIterations;
  }
  return resolved;
}

double boundingBoxDiagonal(const SurfaceMeshData& meshData) {
  if (meshData.vertices.empty()) {
    return 1.0;
  }

  Vec3 minCorner = meshData.vertices.front();
  Vec3 maxCorner = meshData.vertices.front();
  for (const Vec3& vertex : meshData.vertices) {
    minCorner = minCorner.cwiseMin(vertex);
    maxCorner = maxCorner.cwiseMax(vertex);
  }
  const double diagonal = (maxCorner - minCorner).norm();
  return diagonal > 0.0 ? diagonal : 1.0;
}

TraceSettings makeTraceDefaults(const SurfaceMeshData& meshData,
                                const gcs::SurfaceMesh& activeMesh) {
  TraceSettings settings;
  settings.traceLength = 100.0 * boundingBoxDiagonal(meshData);
  settings.maxIterations = std::max<size_t>(100, 100 * activeMesh.nFaces());
  return settings;
}

} // namespace

ReferenceGeometry::ReferenceGeometry(SurfaceMeshData meshData)
    : meshData_(std::move(meshData)), surface_(makeGeometryCentralSurface(meshData_)) {}

SurfaceMeshData& ReferenceGeometry::meshData() {
  return meshData_;
}

const SurfaceMeshData& ReferenceGeometry::meshData() const {
  return meshData_;
}

GeometryCentralSurface& ReferenceGeometry::surface() {
  return surface_;
}

const GeometryCentralSurface& ReferenceGeometry::surface() const {
  return surface_;
}

ActiveIntrinsicDomain::ActiveIntrinsicDomain(ReferenceGeometry& reference,
                                             const IntrinsicConstructionOptions&,
                                             const RefinementOptions& refinementOptions) {
  triangulation_ = std::make_unique<gcs::SignpostIntrinsicTriangulation>(
      *reference.surface().mesh,
      *reference.surface().geometry);

  if (refinementOptions.mode == RefinementMode::DelaunayFlip) {
    triangulation_->flipToDelaunay();
  } else if (refinementOptions.mode == RefinementMode::DelaunayRefine) {
    const double angleThreshold = refinementOptions.angleThreshold.value_or(25.0);
    const double circumradiusThreshold =
        refinementOptions.circumradiusThreshold.value_or(std::numeric_limits<double>::infinity());
    const size_t maxInsertions =
        refinementOptions.maxInsertions.value_or(geometrycentral::INVALID_IND);
    triangulation_->delaunayRefine(angleThreshold, circumradiusThreshold, maxInsertions);
  }
}

gcs::ManifoldSurfaceMesh& ActiveIntrinsicDomain::mesh() {
  return *triangulation_->intrinsicMesh;
}

gcs::IntrinsicGeometryInterface& ActiveIntrinsicDomain::geometry() {
  return *triangulation_;
}

gcs::IntrinsicTriangulation& ActiveIntrinsicDomain::triangulation() {
  return *triangulation_;
}

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
          heatOptions_.diffusionTimeCoefficient)) {}

DrapeResult GeoDrapeSolver::solve(const Vec2& seedXY,
                                  double angleDegrees,
                                  const DrapeSolveOptions& solveOptions) {
  return solve(seedXY, angleDegrees, 90.0, solveOptions);
}

DrapeResult GeoDrapeSolver::solve(const Vec2& seedXY,
                                  double fabricAngleDegrees,
                                  double fiberAngleDegrees,
                                  const DrapeSolveOptions& solveOptions) {
  const TraceSettings trace = resolveTraceSettings(traceDefaults_, solveOptions.advanced.trace);
  const IntrinsicSolveInput input = adaptExtrinsicInput(
      seedXY,
      fabricAngleDegrees,
      fiberAngleDegrees,
      solveOptions.mode,
      trace);
  lastIntrinsicResult_ = solveCore(input);
  DrapeResult result = retrieveFromCore(
      *lastIntrinsicResult_,
      solveOptions.retrieval,
      solveOptions.sampleVertexShear || solveOptions.sampleSecondaryShear);
  if (solveOptions.sampleSecondaryShear && solveOptions.mode == DrapeSolveMode::Complete) {
    result.faceShearAnglesDegrees = lastIntrinsicResult_->faceShearAnglesDegrees;
  }
  return result;
}

DrapeResult GeoDrapeSolver::retrieve(ResultDomain retrieval, bool sampleVertexShear) {
  if (!lastIntrinsicResult_) {
    throw std::runtime_error("GeoDrapeSolver::retrieve() requires a previous solve");
  }
  return retrieveFromCore(*lastIntrinsicResult_, retrieval, sampleVertexShear);
}

IntrinsicSolveInput GeoDrapeSolver::adaptExtrinsicInput(const Vec2& seedXY,
                                                        double fabricAngleDegrees,
                                                        double fiberAngleDegrees,
                                                        DrapeSolveMode mode,
                                                        const TraceSettings& trace) {
  const std::optional<SeedProjection> seed = projectPointXYToMesh(reference_.meshData(), seedXY);
  if (!seed) {
    throw std::runtime_error("GeoDrapeSolver failed to project seed point to mesh");
  }

  IntrinsicSolveInput input;
  input.seed = seed->surfacePoint;
  input.directions = generateFamilyDirections(fabricAngleDegrees, fiberAngleDegrees);
  input.mode = mode;
  input.trace = trace;
  return input;
}

CoreIntrinsicResult GeoDrapeSolver::solveCore(const IntrinsicSolveInput& input) {
  CoreIntrinsicResult result;
  result.mode = input.mode;
  result.directions = input.directions;
  result.seed.surfacePoint = input.seed;

  const Face& face = reference_.meshData().faces.at(input.seed.faceIndex);
  result.seed.cartesian =
      input.seed.barycentric(0) * reference_.meshData().vertices[face[0]] +
      input.seed.barycentric(1) * reference_.meshData().vertices[face[1]] +
      input.seed.barycentric(2) * reference_.meshData().vertices[face[2]];

  result.generators = traceGenerators(
      reference_.surface(),
      result.seed.surfacePoint,
      result.directions,
      input.trace);
  result.sourceCurves = pairOppositeGeneratorTraces(result.generators);

  const bool computeDistances = input.mode != DrapeSolveMode::Fast;
  result.customHeatSolves = customHeatSolver_->solve(result.sourceCurves, heatOptions_, computeDistances);
  if (computeDistances) {
    result.distances = std::array<std::vector<double>, 2>{};
  }

  for (size_t i = 0; i < result.customHeatSolves.size(); ++i) {
    result.faceDirections[i] = result.customHeatSolves[i].normalizedFaceDirections;
    if (computeDistances) {
      (*result.distances)[i] = result.customHeatSolves[i].distance;
    }
  }

  if (input.mode == DrapeSolveMode::Complete) {
    result.gradients = std::array<std::vector<Vec3>, 2>{};
    (*result.gradients)[0] = computeVertexScalarGradients(reference_.meshData(), (*result.distances)[0]);
    (*result.gradients)[1] = computeVertexScalarGradients(reference_.meshData(), (*result.distances)[1]);
    result.faceShearAnglesDegrees = averageVertexScalarsToFaces(
        reference_.meshData(),
        computeShearAnglesDegrees((*result.gradients)[0], (*result.gradients)[1]));
  } else {
    result.faceShearAnglesDegrees = computeFaceShearAnglesDegrees(
        reference_.surface(),
        result.faceDirections[0],
        result.faceDirections[1]);
  }

  return result;
}

DrapeResult GeoDrapeSolver::retrieveFromCore(const CoreIntrinsicResult& core,
                                             ResultDomain retrieval,
                                             bool sampleVertexShear) const {
  if (retrieval != ResultDomain::Extrinsic) {
    throw std::runtime_error("only extrinsic retrieval is implemented in this architecture checkpoint");
  }

  DrapeResult result;
  result.mode = core.mode;
  result.seed = core.seed;
  result.directions = core.directions;
  result.generators = core.generators;
  result.sourceCurves = core.sourceCurves;
  result.customHeatSolves = core.customHeatSolves;
  result.faceDirections = core.faceDirections;
  result.distances = core.distances;
  result.gradients = core.gradients;

  if (core.mode == DrapeSolveMode::Complete) {
    if (core.gradients) {
      result.vertexShearAnglesDegrees =
          computeShearAnglesDegrees((*core.gradients)[0], (*core.gradients)[1]);
    }
  } else {
    result.faceShearAnglesDegrees = core.faceShearAnglesDegrees;
    if (sampleVertexShear && core.faceShearAnglesDegrees) {
      result.vertexShearAnglesDegrees = averageFaceScalarsToVertices(
          reference_.meshData(),
          *core.faceShearAnglesDegrees,
          FaceScalarAveraging::FaceArea);
    }
  }

  return result;
}

DrapeResult solveDrape(const SurfaceMeshData& mesh,
                       const Vec2& seedXY,
                       double angleDegrees,
                       const SignedHeatSolveOptions& heatOptions,
                       const DrapeSolveOptions& solveOptions) {
  GeoDrapeSolver solver(mesh, heatOptions);
  return solver.solve(seedXY, angleDegrees, solveOptions);
}

} // namespace geodesic_draping

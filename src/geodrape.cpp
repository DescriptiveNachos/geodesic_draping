#include "geodesic_draping/geodrape.h"

#include "geometrycentral/surface/signpost_intrinsic_triangulation.h"
#include "geometrycentral/surface/trace_geodesic.h"

#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>

namespace geodesic_draping {
namespace gcs = geometrycentral::surface;

namespace {

constexpr double pi = 3.141592653589793238462643383279502884;

geometrycentral::Vector3 toGcVector3(const Vec3& value) {
  return geometrycentral::Vector3{value.x(), value.y(), value.z()};
}

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
  (void)meshData;
  (void)activeMesh;
  TraceSettings settings;
  return settings;
}

gcs::SurfacePoint toFaceSurfacePoint(gcs::SurfaceMesh& mesh, const BarycentricPoint& point) {
  return gcs::SurfacePoint(
      mesh.face(point.faceIndex),
      geometrycentral::Vector3{
          point.barycentric.x(),
          point.barycentric.y(),
          point.barycentric.z(),
      });
}

gcs::BarycentricVector normalizeVector(gcs::BarycentricVector vector,
                                       gcs::IntrinsicGeometryInterface& geometry) {
  const double magnitude = vector.norm(geometry);
  if (magnitude == 0.0) {
    throw std::runtime_error("cannot normalize zero tangent direction");
  }
  return vector / magnitude;
}

gcs::BarycentricVector faceTangentVectorToBarycentric(gcs::Face face,
                                                      const geometrycentral::Vector2& tangentVector,
                                                      gcs::IntrinsicGeometryInterface& geometry) {
  geometry.requireHalfedgeVectorsInFace();
  const geometrycentral::Vector2 p1 = geometry.halfedgeVectorsInFace[face.halfedge()];
  const geometrycentral::Vector2 p2 = -geometry.halfedgeVectorsInFace[face.halfedge().next().next()];
  const double det = p1.x * p2.y - p1.y * p2.x;
  if (det == 0.0) {
    geometry.unrequireHalfedgeVectorsInFace();
    throw std::runtime_error("cannot convert tangent vector on degenerate face");
  }
  const double b1 = (tangentVector.x * p2.y - tangentVector.y * p2.x) / det;
  const double b2 = (p1.x * tangentVector.y - p1.y * tangentVector.x) / det;
  geometry.unrequireHalfedgeVectorsInFace();
  return gcs::BarycentricVector(face, geometrycentral::Vector3{-b1 - b2, b1, b2});
}

gcs::BarycentricVector embeddedActiveFaceDirection(ReferenceGeometry& reference,
                                                  ActiveIntrinsicDomain& activeDomain,
                                                  gcs::Face activeFace,
                                                  double angleDegrees) {
  std::array<Vec3, 3> positions;
  size_t i = 0;
  for (gcs::Vertex vertex : activeFace.adjacentVertices()) {
    const gcs::SurfacePoint inputPoint = activeDomain.intrinsicToInput(gcs::SurfacePoint(vertex));
    positions[i] = interpolateSurfacePoint(inputPoint, *reference.surface().geometry);
    ++i;
  }

  const Vec3 ambient(std::cos(angleDegrees * pi / 180.0),
                     std::sin(angleDegrees * pi / 180.0),
                     0.0);
  const Vec3 e1 = positions[1] - positions[0];
  const Vec3 e2 = positions[2] - positions[0];
  Eigen::Matrix2d gram;
  gram << e1.dot(e1), e1.dot(e2),
          e2.dot(e1), e2.dot(e2);
  Eigen::Vector2d rhs;
  rhs << ambient.dot(e1), ambient.dot(e2);
  const Eigen::Vector2d uv = gram.ldlt().solve(rhs);
  return normalizeVector(
      gcs::BarycentricVector(activeFace, geometrycentral::Vector3{-uv.x() - uv.y(), uv.x(), uv.y()}),
      activeDomain.geometry());
}

gcs::BarycentricVector intrinsicDirectionFromInputAngle(ReferenceGeometry& reference,
                                                       ActiveIntrinsicDomain& activeDomain,
                                                       const gcs::SurfacePoint& inputSeed,
                                                       const gcs::SurfacePoint& intrinsicSeed,
                                                       double angleDegrees,
                                                       bool preservesInputConnectivity) {
  auto& inputGeometry = *reference.surface().geometry;
  const gcs::SurfacePoint inputFaceSeed = inputSeed.inSomeFace();
  inputGeometry.requireFaceTangentBasis();

  const Vec3 ambientDirection(std::cos(angleDegrees * pi / 180.0),
                              std::sin(angleDegrees * pi / 180.0),
                              0.0);
  const geometrycentral::Vector3 ambient = toGcVector3(ambientDirection);
  const geometrycentral::Vector3 basisX = inputGeometry.faceTangentBasis[inputFaceSeed.face][0];
  const geometrycentral::Vector3 basisY = inputGeometry.faceTangentBasis[inputFaceSeed.face][1];
  geometrycentral::Vector2 tangentDirection{
      dot(ambient, basisX),
      dot(ambient, basisY),
  };
  inputGeometry.unrequireFaceTangentBasis();

  const double tangentNorm = norm(tangentDirection);
  if (tangentNorm == 0.0) {
    throw std::runtime_error("fabric direction is degenerate in the seed tangent plane");
  }
  tangentDirection /= tangentNorm;

  if (preservesInputConnectivity) {
    const gcs::SurfacePoint intrinsicSeedFace = intrinsicSeed.inSomeFace();
    return normalizeVector(
        faceTangentVectorToBarycentric(intrinsicSeedFace.face, tangentDirection, activeDomain.geometry()),
        activeDomain.geometry());
  }

  const double baseStep = std::max(1e-9, 1e-6 * boundingBoxDiagonal(reference.meshData()));
  gcs::TraceOptions options;
  options.includePath = false;
  options.errorOnProblem = false;
  options.maxIters = 16;

  const gcs::SurfacePoint intrinsicSeedFace = intrinsicSeed.inSomeFace();
  for (size_t attempt = 0; attempt < 8; ++attempt) {
    const double step = baseStep * std::pow(0.25, static_cast<double>(attempt));
    const gcs::TraceGeodesicResult inputTrace =
        gcs::traceGeodesic(inputGeometry, inputSeed, step * tangentDirection, options);
    const gcs::SurfacePoint intrinsicEndpoint =
        activeDomain.inputToIntrinsic(inputTrace.endPoint).inSomeFace();
    try {
      const gcs::SurfacePoint endpointInSeedFace = intrinsicEndpoint.inFace(intrinsicSeedFace.face);
      return normalizeVector(gcs::BarycentricVector(intrinsicSeedFace, endpointInSeedFace),
                             activeDomain.geometry());
    } catch (const std::exception&) {
      // Retry with a shorter local displacement if the mapped endpoint crossed a face boundary.
    }
  }

  return embeddedActiveFaceDirection(reference, activeDomain, intrinsicSeedFace.face, angleDegrees);
}

GeneratorTrace traceActiveGenerator(ReferenceGeometry& reference,
                                    ActiveIntrinsicDomain& activeDomain,
                                    const gcs::SurfacePoint& start,
                                    const gcs::BarycentricVector& direction,
                                    const TraceSettings& settings) {
  gcs::SurfacePoint startFace = start.inSomeFace();
  gcs::BarycentricVector directionInFace = direction.inFace(startFace.face);
  directionInFace = normalizeVector(directionInFace, activeDomain.geometry()) * settings.traceLength;

  gcs::TraceOptions options;
  options.includePath = true;
  options.errorOnProblem = false;
  options.barrierEdges = nullptr;
  options.maxIters = settings.maxIterations;

  const gcs::TraceGeodesicResult traced = gcs::traceGeodesic(
      activeDomain.geometry(),
      startFace.face,
      startFace.faceCoords,
      directionInFace.faceCoords,
      options);
  if (!traced.hasPath) {
    throw std::runtime_error("active traceGeodesic did not return a path");
  }

  GeneratorTrace trace;
  trace.hitBoundary = traced.hitBoundary;
  trace.length = traced.length;
  trace.points.reserve(traced.pathPoints.size());
  trace.surfaceReferences.reserve(traced.pathPoints.size());
  for (const gcs::SurfacePoint& intrinsicPoint : traced.pathPoints) {
    trace.surfaceReferences.push_back(toSurfaceReference(intrinsicPoint));
    const gcs::SurfacePoint inputPoint = activeDomain.intrinsicToInput(intrinsicPoint);
    trace.points.push_back(interpolateSurfacePoint(inputPoint, *reference.surface().geometry));
  }
  return trace;
}

std::array<GeneratorTrace, 4> traceActiveGenerators(ReferenceGeometry& reference,
                                                    ActiveIntrinsicDomain& activeDomain,
                                                    const gcs::SurfacePoint& start,
                                                    const std::array<gcs::BarycentricVector, 4>& directions,
                                                    const TraceSettings& settings) {
  return {
      traceActiveGenerator(reference, activeDomain, start, directions[0], settings),
      traceActiveGenerator(reference, activeDomain, start, directions[1], settings),
      traceActiveGenerator(reference, activeDomain, start, directions[2], settings),
      traceActiveGenerator(reference, activeDomain, start, directions[3], settings),
  };
}

FaceHeatDirectionField computeIntrinsicFaceScalarGradients(gcs::SurfaceMesh& mesh,
                                                           gcs::IntrinsicGeometryInterface& geometry,
                                                           const std::vector<double>& scalarField) {
  if (scalarField.size() != mesh.nVertices()) {
    throw std::runtime_error("computeIntrinsicFaceScalarGradients requires one scalar per active vertex");
  }

  FaceHeatDirectionField gradients(mesh.nFaces(), Vec3::Zero());
  for (gcs::Face face : mesh.faces()) {
    gcs::BarycentricVector gradient(face);
    for (gcs::Halfedge halfedge : face.adjacentHalfedges()) {
      const gcs::BarycentricVector edgeVector(halfedge.next(), face);
      const gcs::BarycentricVector edgePerp = edgeVector.rotate90(geometry);
      gradient += edgePerp * scalarField[halfedge.vertex().getIndex()];
    }
    const double magnitude = gradient.norm(geometry);
    if (magnitude > 0.0) {
      gradient /= magnitude;
    }
    gradients[face.getIndex()] =
        Vec3(gradient.faceCoords.x, gradient.faceCoords.y, gradient.faceCoords.z);
  }
  return gradients;
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

gcs::SurfacePoint ActiveIntrinsicDomain::inputToIntrinsic(const gcs::SurfacePoint& pointOnInput) {
  return triangulation_->equivalentPointOnIntrinsic(pointOnInput);
}

gcs::SurfacePoint ActiveIntrinsicDomain::intrinsicToInput(const gcs::SurfacePoint& pointOnIntrinsic) {
  return triangulation_->equivalentPointOnInput(pointOnIntrinsic);
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
      heatOptions_.diffusionTimeCoefficient)),
      preservesInputConnectivity_(refinementOptions.mode == RefinementMode::None) {}

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
  const gcs::SurfacePoint inputSeed = toFaceSurfacePoint(*reference_.surface().mesh, seed->surfacePoint);
  const gcs::SurfacePoint intrinsicSeed = activeDomain_.inputToIntrinsic(inputSeed).inSomeFace();
  const gcs::BarycentricVector direction0 = intrinsicDirectionFromInputAngle(
      reference_,
      activeDomain_,
      inputSeed,
      intrinsicSeed,
      fabricAngleDegrees,
      preservesInputConnectivity_);
  const gcs::BarycentricVector direction1 = intrinsicDirectionFromInputAngle(
      reference_,
      activeDomain_,
      inputSeed,
      intrinsicSeed,
      fabricAngleDegrees + fiberAngleDegrees,
      preservesInputConnectivity_);

  IntrinsicSolveInput input;
  input.seed = intrinsicSeed;
  input.directions = {direction0, -direction0, direction1, -direction1};
  input.cartesianDirections = generateFamilyDirections(fabricAngleDegrees, fiberAngleDegrees);
  input.mode = mode;
  input.trace = trace;
  return input;
}

CoreIntrinsicResult GeoDrapeSolver::solveCore(const IntrinsicSolveInput& input) {
  CoreIntrinsicResult result;
  result.mode = input.mode;
  result.directions = input.cartesianDirections;
  const gcs::SurfacePoint inputSeed = activeDomain_.intrinsicToInput(input.seed).inSomeFace();
  result.seed.surfacePoint.faceIndex = inputSeed.face.getIndex();
  result.seed.surfacePoint.barycentric =
      Vec3(inputSeed.faceCoords.x, inputSeed.faceCoords.y, inputSeed.faceCoords.z);

  const Face& face = reference_.meshData().faces.at(result.seed.surfacePoint.faceIndex);
  result.seed.cartesian =
      result.seed.surfacePoint.barycentric(0) * reference_.meshData().vertices[face[0]] +
      result.seed.surfacePoint.barycentric(1) * reference_.meshData().vertices[face[1]] +
      result.seed.surfacePoint.barycentric(2) * reference_.meshData().vertices[face[2]];

  if (preservesInputConnectivity_) {
    result.generators = traceGenerators(
        reference_.surface(),
        result.seed.surfacePoint,
        result.directions,
        input.trace);
  } else {
    result.generators = traceActiveGenerators(
        reference_,
        activeDomain_,
        input.seed,
        input.directions,
        input.trace);
  }
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
    result.faceDirections[0] = computeIntrinsicFaceScalarGradients(
        activeDomain_.mesh(),
        activeDomain_.geometry(),
        (*result.distances)[0]);
    result.faceDirections[1] = computeIntrinsicFaceScalarGradients(
        activeDomain_.mesh(),
        activeDomain_.geometry(),
        (*result.distances)[1]);
    if ((*result.distances)[0].size() == reference_.meshData().vertices.size()) {
      result.gradients = std::array<std::vector<Vec3>, 2>{};
      (*result.gradients)[0] = computeVertexScalarGradients(reference_.meshData(), (*result.distances)[0]);
      (*result.gradients)[1] = computeVertexScalarGradients(reference_.meshData(), (*result.distances)[1]);
    }
    result.faceShearAnglesDegrees = computeFaceShearAnglesDegrees(
        activeDomain_.mesh(),
        activeDomain_.geometry(),
        result.faceDirections[0],
        result.faceDirections[1]);
  } else {
    result.faceShearAnglesDegrees = computeFaceShearAnglesDegrees(
        activeDomain_.mesh(),
        activeDomain_.geometry(),
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

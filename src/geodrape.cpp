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

TangentVectorRef toTangentVectorRef(const gcs::BarycentricVector& vector) {
  const gcs::BarycentricVector faceVector = vector.inSomeFace();
  TangentVectorRef ref;
  ref.type = SurfaceReferenceType::Face;
  ref.elementIndex = faceVector.face.getIndex();
  ref.coords = {faceVector.faceCoords.x, faceVector.faceCoords.y, faceVector.faceCoords.z};
  return ref;
}

std::vector<Face> meshFaces(gcs::SurfaceMesh& mesh) {
  std::vector<Face> faces;
  faces.reserve(mesh.nFaces());
  for (gcs::Face face : mesh.faces()) {
    Face out{};
    size_t i = 0;
    for (gcs::Vertex vertex : face.adjacentVertices()) {
      out[i] = vertex.getIndex();
      ++i;
    }
    faces.push_back(out);
  }
  return faces;
}

size_t localSideIndex(gcs::Halfedge target) {
  size_t side = 0;
  for (gcs::Halfedge halfedge : target.face().adjacentHalfedges()) {
    if (halfedge == target) {
      return side;
    }
    ++side;
  }
  throw std::runtime_error("halfedge is not incident on its face");
}

ResultMesh makeExtrinsicResultMesh(const SurfaceMeshData& meshData) {
  ResultMesh mesh;
  mesh.domain = ResultDomain::Extrinsic;
  mesh.faces = meshData.faces;
  mesh.vertices3D = meshData.vertices;
  return mesh;
}

ResultMesh makeIntrinsicResultMesh(gcs::SurfaceMesh& mesh,
                                   gcs::IntrinsicGeometryInterface& geometry) {
  ResultMesh out;
  out.domain = ResultDomain::Intrinsic;
  out.faces = meshFaces(mesh);
  out.edgeLengths = std::vector<std::array<double, 3>>{};
  out.gluingMap = std::vector<FaceGluingMap>{};
  out.edgeLengths->reserve(mesh.nFaces());
  out.gluingMap->reserve(mesh.nFaces());

  geometry.requireEdgeLengths();
  for (gcs::Face face : mesh.faces()) {
    std::array<double, 3> lengths{};
    FaceGluingMap gluing{};
    size_t side = 0;
    for (gcs::Halfedge halfedge : face.adjacentHalfedges()) {
      lengths[side] = geometry.edgeLengths[halfedge.edge()];
      if (halfedge.twin().isInterior()) {
        const gcs::Halfedge twin = halfedge.twin();
        gluing[side] = {
            static_cast<int>(twin.face().getIndex()),
            static_cast<int>(localSideIndex(twin)),
        };
      } else {
        gluing[side] = {-1, -1};
      }
      ++side;
    }
    out.edgeLengths->push_back(lengths);
    out.gluingMap->push_back(gluing);
  }
  geometry.unrequireEdgeLengths();
  return out;
}

std::vector<Vec3> toVec3Vector(const gcs::VertexData<geometrycentral::Vector3>& values) {
  std::vector<Vec3> out(values.getMesh()->nVertices(), Vec3::Zero());
  for (gcs::Vertex vertex : values.getMesh()->vertices()) {
    const geometrycentral::Vector3& value = values[vertex];
    out[vertex.getIndex()] = Vec3(value.x, value.y, value.z);
  }
  return out;
}

std::vector<double> vertexDataToVector(const gcs::VertexData<double>& values) {
  std::vector<double> out(values.getMesh()->nVertices(), 0.0);
  for (gcs::Vertex vertex : values.getMesh()->vertices()) {
    out[vertex.getIndex()] = values[vertex];
  }
  return out;
}

std::vector<double> faceDataToVector(const gcs::FaceData<double>& values) {
  std::vector<double> out(values.getMesh()->nFaces(), 0.0);
  for (gcs::Face face : values.getMesh()->faces()) {
    out[face.getIndex()] = values[face];
  }
  return out;
}

FaceHeatDirectionField faceVectorDataToVector(const gcs::FaceData<Vec3>& values) {
  FaceHeatDirectionField out(values.getMesh()->nFaces(), Vec3::Zero());
  for (gcs::Face face : values.getMesh()->faces()) {
    out[face.getIndex()] = values[face];
  }
  return out;
}

ResultMesh makeSubdivisionResultMesh(gcs::CommonSubdivision& subdivision,
                                     gcs::VertexPositionGeometry& inputGeometry) {
  subdivision.constructMesh();
  ResultMesh out;
  out.domain = ResultDomain::Subdivision;
  out.faces = meshFaces(*subdivision.mesh);
  inputGeometry.requireVertexPositions();
  out.vertices3D = toVec3Vector(subdivision.interpolateAcrossA(inputGeometry.vertexPositions));
  inputGeometry.unrequireVertexPositions();
  return out;
}

DrapeTrace makeTrace(const GeneratorTrace& generator, ResultDomain domain) {
  DrapeTrace trace;
  trace.hitBoundary = generator.hitBoundary;
  trace.length = generator.length;
  if (domain == ResultDomain::Intrinsic) {
    trace.intrinsicPoints = generator.surfaceReferences;
  } else {
    trace.extrinsicPoints = generator.points;
  }
  return trace;
}

std::array<TraceFamily, 2> makeTraceFamilies(const std::array<GeneratorTrace, 4>& generators,
                                             ResultDomain domain) {
  return {
      TraceFamily{makeTrace(generators[0], domain), makeTrace(generators[1], domain)},
      TraceFamily{makeTrace(generators[2], domain), makeTrace(generators[3], domain)},
  };
}

std::vector<double> averageIntrinsicFaceScalarsToVertices(gcs::SurfaceMesh& mesh,
                                                          const std::vector<double>& faceScalars) {
  if (faceScalars.size() != mesh.nFaces()) {
    throw std::runtime_error("averageIntrinsicFaceScalarsToVertices requires one scalar per active face");
  }

  std::vector<double> accumulated(mesh.nVertices(), 0.0);
  std::vector<double> counts(mesh.nVertices(), 0.0);
  for (gcs::Face face : mesh.faces()) {
    for (gcs::Vertex vertex : face.adjacentVertices()) {
      accumulated[vertex.getIndex()] += faceScalars[face.getIndex()];
      counts[vertex.getIndex()] += 1.0;
    }
  }
  for (size_t i = 0; i < accumulated.size(); ++i) {
    if (counts[i] > 0.0) {
      accumulated[i] /= counts[i];
    }
  }
  return accumulated;
}

std::vector<double> toVector(const gcs::VertexData<double>& values) {
  return vertexDataToVector(values);
}

std::vector<double> restrictVertexScalarsToInput(gcs::IntrinsicTriangulation& triangulation,
                                                 const std::vector<double>& valuesOnIntrinsic) {
  if (valuesOnIntrinsic.size() != triangulation.intrinsicMesh->nVertices()) {
    throw std::runtime_error("restrictVertexScalarsToInput requires one scalar per active intrinsic vertex");
  }
  gcs::VertexData<double> activeValues(*triangulation.intrinsicMesh, 0.0);
  for (gcs::Vertex vertex : triangulation.intrinsicMesh->vertices()) {
    activeValues[vertex] = valuesOnIntrinsic[vertex.getIndex()];
  }
  return toVector(triangulation.restrictToInput(activeValues));
}

gcs::VertexData<double> activeVertexData(gcs::SurfaceMesh& mesh, const std::vector<double>& values) {
  if (values.size() != mesh.nVertices()) {
    throw std::runtime_error("activeVertexData requires one scalar per active vertex");
  }
  gcs::VertexData<double> data(mesh, 0.0);
  for (gcs::Vertex vertex : mesh.vertices()) {
    data[vertex] = values[vertex.getIndex()];
  }
  return data;
}

gcs::FaceData<double> activeFaceData(gcs::SurfaceMesh& mesh, const std::vector<double>& values) {
  if (values.size() != mesh.nFaces()) {
    throw std::runtime_error("activeFaceData requires one scalar per active face");
  }
  gcs::FaceData<double> data(mesh, 0.0);
  for (gcs::Face face : mesh.faces()) {
    data[face] = values[face.getIndex()];
  }
  return data;
}

gcs::FaceData<Vec3> activeFaceVectorData(gcs::SurfaceMesh& mesh, const FaceHeatDirectionField& values) {
  if (values.size() != mesh.nFaces()) {
    throw std::runtime_error("activeFaceVectorData requires one vector per active face");
  }
  Eigen::Matrix<Vec3, Eigen::Dynamic, 1> initial(static_cast<Eigen::Index>(mesh.nFaces()));
  for (Eigen::Index i = 0; i < initial.rows(); ++i) {
    initial[i] = Vec3::Zero();
  }
  gcs::FaceData<Vec3> data(mesh, initial);
  for (gcs::Face face : mesh.faces()) {
    data[face] = values[face.getIndex()];
  }
  return data;
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
  result.intrinsicSeed = toSurfaceReference(input.seed);
  result.intrinsicDirections = {
      toTangentVectorRef(input.directions[0]),
      toTangentVectorRef(input.directions[1]),
      toTangentVectorRef(input.directions[2]),
      toTangentVectorRef(input.directions[3]),
  };
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
                                             bool sampleVertexShear) {
  DrapeResult result;
  result.domain = retrieval;
  result.mode = core.mode;

  if (retrieval == ResultDomain::Intrinsic) {
    result.mesh = makeIntrinsicResultMesh(activeDomain_.mesh(), activeDomain_.geometry());
    result.origin.intrinsicPoint = core.intrinsicSeed;
    result.origin.intrinsicFamilyDirections = {core.intrinsicDirections[0], core.intrinsicDirections[2]};
    result.traces = makeTraceFamilies(core.generators, ResultDomain::Intrinsic);
    result.faceDirections = core.faceDirections;
    result.faceShearAnglesDegrees = core.faceShearAnglesDegrees;
    result.distances = core.distances;
    if (sampleVertexShear && core.faceShearAnglesDegrees) {
      result.vertexShearAnglesDegrees =
          averageIntrinsicFaceScalarsToVertices(activeDomain_.mesh(), *core.faceShearAnglesDegrees);
    }
  } else if (retrieval == ResultDomain::Subdivision) {
    gcs::CommonSubdivision& subdivision = activeDomain_.triangulation().getCommonSubdivision();
    subdivision.constructMesh();
    result.mesh = makeSubdivisionResultMesh(subdivision, *reference_.surface().geometry);
    result.origin.intrinsicPoint = core.intrinsicSeed;
    result.origin.intrinsicFamilyDirections = {core.intrinsicDirections[0], core.intrinsicDirections[2]};
    result.origin.extrinsicPoint = core.seed.cartesian;
    result.origin.extrinsicFamilyDirections = {core.directions[0], core.directions[2]};
    result.traces = makeTraceFamilies(core.generators, ResultDomain::Subdivision);

    result.faceDirections = {
        faceVectorDataToVector(subdivision.copyFromB(
            activeFaceVectorData(activeDomain_.mesh(), core.faceDirections[0]))),
        faceVectorDataToVector(subdivision.copyFromB(
            activeFaceVectorData(activeDomain_.mesh(), core.faceDirections[1]))),
    };
    if (core.faceShearAnglesDegrees) {
      result.faceShearAnglesDegrees = faceDataToVector(subdivision.copyFromB(
          activeFaceData(activeDomain_.mesh(), *core.faceShearAnglesDegrees)));
    }
    if (core.distances) {
      result.distances = std::array<std::vector<double>, 2>{
          vertexDataToVector(subdivision.interpolateAcrossB(
              activeVertexData(activeDomain_.mesh(), (*core.distances)[0]))),
          vertexDataToVector(subdivision.interpolateAcrossB(
              activeVertexData(activeDomain_.mesh(), (*core.distances)[1]))),
      };
    }
    if (sampleVertexShear && core.faceShearAnglesDegrees) {
      result.vertexShearAnglesDegrees =
          averageFaceScalarsToVertices(
              SurfaceMeshData{*result.mesh.vertices3D, result.mesh.faces},
              *result.faceShearAnglesDegrees,
              FaceScalarAveraging::FaceArea);
    }
  } else {
    result.mesh = makeExtrinsicResultMesh(reference_.meshData());
    result.origin.intrinsicPoint = core.intrinsicSeed;
    result.origin.intrinsicFamilyDirections = {core.intrinsicDirections[0], core.intrinsicDirections[2]};
    result.origin.extrinsicPoint = core.seed.cartesian;
    result.origin.extrinsicFamilyDirections = {core.directions[0], core.directions[2]};
    result.traces = makeTraceFamilies(core.generators, ResultDomain::Extrinsic);
    result.faceDirections = core.faceDirections;

    if (core.distances) {
      if ((*core.distances)[0].size() == reference_.meshData().vertices.size()) {
        result.distances = core.distances;
      } else {
        result.distances = std::array<std::vector<double>, 2>{
            restrictVertexScalarsToInput(activeDomain_.triangulation(), (*core.distances)[0]),
            restrictVertexScalarsToInput(activeDomain_.triangulation(), (*core.distances)[1]),
        };
      }
    }

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

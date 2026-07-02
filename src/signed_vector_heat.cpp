#include "geodesic_draping/signed_vector_heat.h"

#include "geometrycentral/numerical/linear_algebra_utilities.h"
#include "geometrycentral/surface/barycentric_vector.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <stdexcept>

namespace geodesic_draping {
namespace gc = geometrycentral;
namespace gcs = geometrycentral::surface;

namespace {

gc::Vector<std::complex<double>> toEigenVector(const EdgeHeatField& values) {
  gc::Vector<std::complex<double>> out(values.size());
  for (size_t i = 0; i < values.size(); ++i) {
    out[static_cast<Eigen::Index>(i)] = values[i];
  }
  return out;
}

EdgeHeatField toStdVector(const gc::Vector<std::complex<double>>& values) {
  EdgeHeatField out(static_cast<size_t>(values.rows()));
  for (Eigen::Index i = 0; i < values.rows(); ++i) {
    out[static_cast<size_t>(i)] = values[i];
  }
  return out;
}

std::vector<SurfaceReference> toSurfaceReferences(const gcs::Curve& curve) {
  std::vector<SurfaceReference> refs;
  refs.reserve(curve.nodes.size());
  for (const gcs::SurfacePoint& point : curve.nodes) {
    refs.push_back(toSurfaceReference(point));
  }
  return refs;
}

Vec3 fromGcVector3(const gc::Vector3& value) {
  return Vec3(value.x, value.y, value.z);
}

Eigen::Index eigenIndex(size_t value) {
  return static_cast<Eigen::Index>(value);
}

double cornerAngle(const gc::Vector3& center, const gc::Vector3& pA, const gc::Vector3& pB) {
  const gc::Vector3 u = pA - center;
  const gc::Vector3 v = pB - center;
  const double denom = norm(u) * norm(v);
  if (denom == 0.0) {
    return 0.0;
  }
  return std::acos(std::clamp(dot(u, v) / denom, -1.0, 1.0));
}

gc::Vector3 faceVectorToExtrinsic(const gcs::Face& face,
                                  const Vec3& coords,
                                  gcs::VertexPositionGeometry& geometry) {
  gc::Vector3 vector{0.0, 0.0, 0.0};
  size_t localIndex = 0;
  for (gcs::Vertex vertex : face.adjacentVertices()) {
    const double coefficient = localIndex == 0 ? coords.x() : (localIndex == 1 ? coords.y() : coords.z());
    vector += coefficient * geometry.vertexPositions[vertex];
    ++localIndex;
  }
  return vector;
}

gcs::SurfaceMesh& requireMesh(GeometryCentralSurface& surface) {
  if (!surface.mesh) {
    throw std::runtime_error("CustomSignedHeatSolver requires a valid geometry-central mesh");
  }
  return *surface.mesh;
}

gcs::IntrinsicGeometryInterface& requireGeometry(GeometryCentralSurface& surface) {
  if (!surface.geometry) {
    throw std::runtime_error("CustomSignedHeatSolver requires a valid geometry-central geometry");
  }
  return *surface.geometry;
}

} // namespace

CustomSignedHeatSolver::CustomSignedHeatSolver(GeometryCentralSurface& surface,
                                               double diffusionTimeCoefficient)
    : mesh_(requireMesh(surface)), geometry_(requireGeometry(surface)) {
  geometry_.requireEdgeLengths();
  geometry_.requireCrouzeixRaviartMassMatrix();

  double meanEdgeLength = 0.0;
  for (gcs::Edge edge : mesh_.edges()) {
    meanEdgeLength += geometry_.edgeLengths[edge];
  }
  meanEdgeLength /= static_cast<double>(mesh_.nEdges());
  meanNodeDistance_ = 0.5 * meanEdgeLength;
  shortTime_ = diffusionTimeCoefficient * meanNodeDistance_ * meanNodeDistance_;

  massMatrix_ = geometry_.crouzeixRaviartMassMatrix;
  doubleMassMatrix_ = buildCrouzeixRaviartDoubleMassMatrix();
  doubleConnectionLaplacian_ = buildCrouzeixRaviartDoubleConnectionLaplacian();
  doubleVectorOperator_ = doubleMassMatrix_ + shortTime_ * doubleConnectionLaplacian_;

  geometry_.unrequireCrouzeixRaviartMassMatrix();
  geometry_.unrequireEdgeLengths();
}

DiffusedHeatFieldResult CustomSignedHeatSolver::solveDiffusedEdgeHeatField(
    const std::vector<SurfaceReference>& sourceCurve,
    const SignedHeatSolveOptions& options) {
  gcs::Curve curve = toGeometryCentralCurve(mesh_, sourceCurve, true);
  const std::vector<gcs::Curve> preprocessedCurves = preprocessCurves({curve});

  DiffusedHeatFieldResult result;
  result.preprocessedSourceCurves.reserve(preprocessedCurves.size());
  for (const gcs::Curve& preprocessedCurve : preprocessedCurves) {
    result.preprocessedSourceCurves.push_back(toSurfaceReferences(preprocessedCurve));
  }
  result.sourceEdgeHeatField = buildSourceEdgeHeatField(preprocessedCurves);
  result.diffusedEdgeHeatField =
      diffuseEdgeHeatField(result.sourceEdgeHeatField, preprocessedCurves, options);
  return result;
}

std::vector<gcs::Curve> CustomSignedHeatSolver::preprocessCurves(
    const std::vector<gcs::Curve>& curves) const {
  std::vector<gcs::Curve> newCurves;
  for (const gcs::Curve& curve : curves) {
    if (curve.nodes.empty()) {
      continue;
    }
    newCurves.emplace_back();
    newCurves.back().isSigned = curve.isSigned;
    const size_t nNodes = curve.nodes.size();
    for (size_t i = 0; i + 1 < nNodes; ++i) {
      const gcs::SurfacePoint& pA = curve.nodes[i];
      const gcs::SurfacePoint& pB = curve.nodes[i + 1];
      newCurves.back().nodes.push_back(pA);
      if (sharedFace(pA, pB) == gcs::Face()) {
        newCurves.emplace_back();
        newCurves.back().isSigned = curve.isSigned;
      }
    }
    newCurves.back().nodes.push_back(curve.nodes[nNodes - 1]);
  }
  return newCurves;
}

EdgeHeatField CustomSignedHeatSolver::buildSourceEdgeHeatField(
    const std::vector<gcs::Curve>& curves) const {
  geometry_.requireEdgeIndices();
  EdgeHeatField source(mesh_.nEdges(), std::complex<double>(0.0, 0.0));
  for (const gcs::Curve& curve : curves) {
    if (!curve.isSigned) {
      throw std::runtime_error("unsigned source curves are not implemented in CustomSignedHeatSolver");
    }
    buildSignedCurveSource(curve, source);
  }
  geometry_.unrequireEdgeIndices();
  return source;
}

EdgeHeatField CustomSignedHeatSolver::diffuseEdgeHeatField(
    const EdgeHeatField& sourceEdgeHeatField,
    const std::vector<gcs::Curve>& curves,
    const SignedHeatSolveOptions& options) {
  ensureHaveHeatFieldSolver();

  const size_t E = mesh_.nEdges();
  const gc::Vector<std::complex<double>> source = toEigenVector(sourceEdgeHeatField);
  gc::Vector<std::complex<double>> diffused;

  if (!options.preserveSourceNormals) {
    diffused = heatFieldSolver_->solve(source);
    return toStdVector(diffused);
  }

  geometry_.requireEdgeIndices();
  gc::Vector<double> rhs = gc::Vector<double>::Zero(2 * E);
  for (size_t i = 0; i < E; ++i) {
    rhs[static_cast<Eigen::Index>(i)] = std::real(source[static_cast<Eigen::Index>(i)]);
    rhs[static_cast<Eigen::Index>(E + i)] = std::imag(source[static_cast<Eigen::Index>(i)]);
  }

  std::vector<Eigen::Triplet<double>> triplets;
  size_t m = 0;
  for (const gcs::Curve& curve : curves) {
    for (size_t i = 0; i + 1 < curve.nodes.size(); ++i) {
      const gcs::SurfacePoint& pA = curve.nodes[i];
      const gcs::SurfacePoint& pB = curve.nodes[i + 1];
      const gcs::Edge commonEdge = sharedEdge(pA, pB);
      if (commonEdge != gcs::Edge()) {
        triplets.emplace_back(static_cast<Eigen::Index>(m),
                              static_cast<Eigen::Index>(geometry_.edgeIndices[commonEdge]),
                              1.0);
        ++m;
        continue;
      }

      const gcs::Face face = sharedFace(pA, pB);
      if (face == gcs::Face()) {
        throw std::runtime_error("source curve segment is not contained in a face");
      }
      const gcs::SurfacePoint midpoint = midSegmentSurfacePoint(pA, pB);
      gcs::BarycentricVector segmentTangent(pA, pB);
      for (gcs::Edge edge : face.adjacentEdges()) {
        const size_t eIdx = geometry_.edgeIndices[edge];
        const double weight = scalarCrouzeixRaviart(midpoint, edge);
        gcs::BarycentricVector edgeVector(edge.halfedge(), face);
        edgeVector /= edgeVector.norm(geometry_);
        const gcs::BarycentricVector edgeNormal = edgeVector.rotate90(geometry_);
        triplets.emplace_back(static_cast<Eigen::Index>(m),
                              static_cast<Eigen::Index>(eIdx),
                              weight * dot(geometry_, edgeVector, segmentTangent));
        triplets.emplace_back(static_cast<Eigen::Index>(m),
                              static_cast<Eigen::Index>(E + eIdx),
                              weight * dot(geometry_, edgeNormal, segmentTangent));
      }
      ++m;
    }
  }

  gc::SparseMatrix<double> constraints(static_cast<Eigen::Index>(m), static_cast<Eigen::Index>(2 * E));
  constraints.setFromTriplets(triplets.begin(), triplets.end());
  gc::SparseMatrix<double> zero(static_cast<Eigen::Index>(m), static_cast<Eigen::Index>(m));
  gc::SparseMatrix<double> lhsTop = gc::horizontalStack<double>({doubleVectorOperator_, constraints.transpose()});
  gc::SparseMatrix<double> lhsBottom = gc::horizontalStack<double>({constraints, zero});
  gc::SparseMatrix<double> lhs = gc::verticalStack<double>({lhsTop, lhsBottom});
  gc::Vector<double> fullRhs = gc::Vector<double>::Zero(static_cast<Eigen::Index>(2 * E + m));
  fullRhs.head(static_cast<Eigen::Index>(2 * E)) = rhs;
  const gc::Vector<double> solution = gc::solveSquare(lhs, fullRhs);
  const gc::Vector<double> x = solution.head(static_cast<Eigen::Index>(2 * E));

  diffused = gc::Vector<std::complex<double>>::Zero(static_cast<Eigen::Index>(E));
  for (size_t i = 0; i < E; ++i) {
    diffused[static_cast<Eigen::Index>(i)] =
        std::complex<double>(x[static_cast<Eigen::Index>(i)], x[static_cast<Eigen::Index>(E + i)]);
  }
  geometry_.unrequireEdgeIndices();
  return toStdVector(diffused);
}

void CustomSignedHeatSolver::ensureHaveHeatFieldSolver() {
  if (heatFieldSolver_) {
    return;
  }

  geometry_.requireCrouzeixRaviartConnectionLaplacian();
  gc::SparseMatrix<std::complex<double>>& connectionLaplacian =
      geometry_.crouzeixRaviartConnectionLaplacian;
  gc::SparseMatrix<std::complex<double>> vectorOperator =
      massMatrix_.cast<std::complex<double>>() + shortTime_ * connectionLaplacian;

  geometry_.requireEdgeCotanWeights();
  bool isDelaunay = true;
  for (gcs::Edge edge : mesh_.edges()) {
    if (geometry_.edgeCotanWeights[edge] < -1e-6) {
      isDelaunay = false;
      break;
    }
  }
  geometry_.unrequireEdgeCotanWeights();

  if (isDelaunay) {
    heatFieldSolver_.reset(new gc::PositiveDefiniteSolver<std::complex<double>>(vectorOperator));
  } else {
    heatFieldSolver_.reset(new gc::SquareSolver<std::complex<double>>(vectorOperator));
  }
}

gc::SparseMatrix<double> CustomSignedHeatSolver::buildCrouzeixRaviartDoubleConnectionLaplacian() const {
  geometry_.requireEdgeIndices();
  geometry_.requireEdgeLengths();
  geometry_.requireHalfedgeCotanWeights();
  geometry_.requireFaceAreas();

  const size_t E = mesh_.nEdges();
  gc::SparseMatrix<double> laplacian(static_cast<Eigen::Index>(2 * E), static_cast<Eigen::Index>(2 * E));
  std::vector<Eigen::Triplet<double>> triplets;
  for (gcs::Face face : mesh_.faces()) {
    for (gcs::Halfedge halfedge : face.adjacentHalfedges()) {
      const gcs::Halfedge heA = halfedge.next();
      const gcs::Halfedge heB = heA.next();
      const size_t iE_i = geometry_.edgeIndices[heA.edge()];
      const size_t iE_j = geometry_.edgeIndices[heB.edge()];

      const double weight = 4.0 * geometry_.halfedgeCotanWeights[halfedge];
      const double sign = (heA.orientation() == heB.orientation()) ? 1.0 : -1.0;
      const double lOpp = geometry_.edgeLengths[halfedge.edge()];
      const double lA = geometry_.edgeLengths[heB.edge()];
      const double lB = geometry_.edgeLengths[heA.edge()];
      const double cosTheta = (lA * lA + lB * lB - lOpp * lOpp) / (2.0 * lA * lB);
      const double sinTheta = 2.0 * geometry_.faceAreas[face] / (lA * lB);
      const double a = weight * sign * cosTheta;
      const double b = weight * sign * sinTheta;

      triplets.emplace_back(eigenIndex(iE_i), eigenIndex(iE_i), weight);
      triplets.emplace_back(eigenIndex(iE_j), eigenIndex(iE_j), weight);
      triplets.emplace_back(eigenIndex(E + iE_i), eigenIndex(E + iE_i), weight);
      triplets.emplace_back(eigenIndex(E + iE_j), eigenIndex(E + iE_j), weight);
      triplets.emplace_back(eigenIndex(iE_i), eigenIndex(iE_j), a);
      triplets.emplace_back(eigenIndex(iE_j), eigenIndex(iE_i), a);
      triplets.emplace_back(eigenIndex(E + iE_i), eigenIndex(E + iE_j), a);
      triplets.emplace_back(eigenIndex(E + iE_j), eigenIndex(E + iE_i), a);
      triplets.emplace_back(eigenIndex(iE_i), eigenIndex(E + iE_j), b);
      triplets.emplace_back(eigenIndex(E + iE_j), eigenIndex(iE_i), b);
      triplets.emplace_back(eigenIndex(iE_j), eigenIndex(E + iE_i), -b);
      triplets.emplace_back(eigenIndex(E + iE_i), eigenIndex(iE_j), -b);
    }
  }
  laplacian.setFromTriplets(triplets.begin(), triplets.end());

  geometry_.unrequireFaceAreas();
  geometry_.unrequireEdgeIndices();
  geometry_.unrequireEdgeLengths();
  geometry_.unrequireHalfedgeCotanWeights();
  return laplacian;
}

gc::SparseMatrix<double> CustomSignedHeatSolver::buildCrouzeixRaviartDoubleMassMatrix() const {
  geometry_.requireEdgeIndices();
  geometry_.requireFaceAreas();

  const size_t E = mesh_.nEdges();
  gc::SparseMatrix<double> mass(static_cast<Eigen::Index>(2 * E), static_cast<Eigen::Index>(2 * E));
  std::vector<Eigen::Triplet<double>> triplets;
  for (gcs::Edge edge : mesh_.edges()) {
    const size_t eIdx = geometry_.edgeIndices[edge];
    double area = 0.0;
    for (gcs::Face face : edge.adjacentFaces()) {
      area += geometry_.faceAreas[face];
    }
    const double value = area / 3.0;
    triplets.emplace_back(eigenIndex(eIdx), eigenIndex(eIdx), value);
    triplets.emplace_back(eigenIndex(E + eIdx), eigenIndex(E + eIdx), value);
  }
  mass.setFromTriplets(triplets.begin(), triplets.end());

  geometry_.unrequireEdgeIndices();
  geometry_.unrequireFaceAreas();
  return mass;
}

void CustomSignedHeatSolver::buildSignedCurveSource(const gcs::Curve& curve,
                                                    EdgeHeatField& sourceEdgeHeatField) const {
  for (size_t i = 0; i + 1 < curve.nodes.size(); ++i) {
    const gcs::SurfacePoint& pA = curve.nodes[i];
    const gcs::SurfacePoint& pB = curve.nodes[i + 1];
    const gcs::Edge commonEdge = sharedEdge(pA, pB);
    if (commonEdge != gcs::Edge()) {
      const size_t eIdx = geometry_.edgeIndices[commonEdge];
      std::complex<double> innerProduct(0.0, 1.0);
      if (pA.vertex == commonEdge.secondVertex()) {
        innerProduct.imag(-1.0);
      }
      sourceEdgeHeatField[eIdx] += lengthOfSegment(pA, pB) * innerProduct;
      continue;
    }

    const gcs::Face commonFace = sharedFace(pA, pB);
    if (commonFace == gcs::Face()) {
      throw std::runtime_error("source curve segment is not contained in a face");
    }
    for (gcs::Edge edge : commonFace.adjacentEdges()) {
      sourceEdgeHeatField[geometry_.edgeIndices[edge]] += projectedNormal(pA, pB, edge);
    }
  }
}

double CustomSignedHeatSolver::lengthOfSegment(const gcs::SurfacePoint& pA,
                                               const gcs::SurfacePoint& pB) const {
  gcs::BarycentricVector segment(pA, pB);
  return segment.norm(geometry_);
}

gcs::SurfacePoint CustomSignedHeatSolver::midSegmentSurfacePoint(const gcs::SurfacePoint& pA,
                                                                 const gcs::SurfacePoint& pB) const {
  const gcs::Face commonFace = sharedFace(pA, pB);
  const gcs::SurfacePoint pAInFace = pA.inFace(commonFace);
  const gcs::SurfacePoint pBInFace = pB.inFace(commonFace);
  return gcs::SurfacePoint(commonFace, 0.5 * (pAInFace.faceCoords + pBInFace.faceCoords));
}

std::complex<double> CustomSignedHeatSolver::projectedNormal(const gcs::SurfacePoint& pA,
                                                             const gcs::SurfacePoint& pB,
                                                             const gcs::Edge& edge) const {
  gcs::BarycentricVector segment(pA, pB);
  const gcs::BarycentricVector segmentNormal = segment.rotate90(geometry_);
  const gcs::Face face = segment.face;

  gc::Vector3 edgeVectorCoords = {0.0, 0.0, 0.0};
  int vertexIndex = 0;
  bool orientation = true;
  for (gcs::Halfedge halfedge : face.adjacentHalfedges()) {
    if (halfedge.edge() == edge) {
      edgeVectorCoords[vertexIndex] = -1.0;
      edgeVectorCoords[(vertexIndex + 1) % 3] = 1.0;
      orientation = (halfedge.vertex() == edge.firstVertex());
      break;
    }
    ++vertexIndex;
  }
  if (!orientation) {
    edgeVectorCoords *= -1.0;
  }
  gcs::BarycentricVector edgeVector(face, edgeVectorCoords);
  edgeVector /= geometry_.edgeLengths[edge];

  const double sinTheta = dot(geometry_, segment, edgeVector);
  const double cosTheta = dot(geometry_, segmentNormal, edgeVector);
  return {cosTheta, sinTheta};
}

double CustomSignedHeatSolver::scalarCrouzeixRaviart(const gcs::SurfacePoint& point,
                                                     const gcs::Edge& edge) const {
  gcs::Face shared = gcs::Face();
  switch (point.type) {
  case gcs::SurfacePointType::Vertex:
    for (gcs::Face face : edge.adjacentFaces()) {
      for (gcs::Vertex vertex : face.adjacentVertices()) {
        if (vertex == point.vertex) {
          shared = face;
          break;
        }
      }
    }
    if (shared == gcs::Face()) return 0.0;
    if (point.vertex == edge.firstVertex() || point.vertex == edge.secondVertex()) return 1.0;
    return -1.0;
  case gcs::SurfacePointType::Edge:
    if (point.edge == edge) return 1.0;
    for (gcs::Face face : edge.adjacentFaces()) {
      for (gcs::Edge adjacentEdge : face.adjacentEdges()) {
        if (adjacentEdge == point.edge) {
          shared = face;
          break;
        }
      }
    }
    if (shared == gcs::Face()) return 0.0;
    {
      const gcs::SurfacePoint pointInFace = point.inFace(shared);
      int index = 0;
      for (gcs::Halfedge halfedge : shared.adjacentHalfedges()) {
        if (halfedge.next().edge() == edge) {
          return 1.0 - 2.0 * pointInFace.faceCoords[index];
        }
        ++index;
      }
    }
    break;
  case gcs::SurfacePointType::Face:
    for (gcs::Face face : edge.adjacentFaces()) {
      if (face == point.face) {
        shared = face;
        break;
      }
    }
    if (shared == gcs::Face()) return 0.0;
    {
      size_t index = 0;
      double value = 0.0;
      for (gcs::Halfedge halfedge : shared.adjacentHalfedges()) {
        const double vertexValue =
            (halfedge.vertex() == edge.firstVertex() || halfedge.vertex() == edge.secondVertex()) ? 1.0 : -1.0;
        value += vertexValue * point.faceCoords[index];
        ++index;
      }
      return value;
    }
  }
  throw std::runtime_error("scalarCrouzeixRaviart received an unsupported SurfacePoint type");
}

FaceHeatDirectionField sampleAndNormalizeFaceDirections(GeometryCentralSurface& surface,
                                                        const EdgeHeatField& diffusedEdgeHeatField) {
  auto& mesh = *surface.mesh;
  auto& geometry = *surface.geometry;
  if (diffusedEdgeHeatField.size() != mesh.nEdges()) {
    throw std::runtime_error("sampleAndNormalizeFaceDirections requires one heat-field value per edge");
  }

  geometry.requireEdgeIndices();
  FaceHeatDirectionField normalized(mesh.nFaces(), Vec3::Zero());
  for (gcs::Face face : mesh.faces()) {
    gc::Vector3 faceCoords = {0.0, 0.0, 0.0};
    for (gcs::Halfedge halfedge : face.adjacentHalfedges()) {
      const size_t eIdx = geometry.edgeIndices[halfedge.edge()];
      gcs::BarycentricVector e1(halfedge, face);
      if (!halfedge.orientation()) {
        e1 *= -1.0;
      }
      gcs::BarycentricVector e2 = e1.rotate90(geometry);
      e1 /= e1.norm(geometry);
      e2 /= e2.norm(geometry);
      faceCoords += std::real(diffusedEdgeHeatField[eIdx]) * e1.faceCoords;
      faceCoords += std::imag(diffusedEdgeHeatField[eIdx]) * e2.faceCoords;
    }
    gcs::BarycentricVector vector(face, faceCoords);
    const double magnitude = vector.norm(geometry);
    if (magnitude != 0.0) {
      faceCoords = (vector / magnitude).faceCoords;
    } else {
      faceCoords = {0.0, 0.0, 0.0};
    }
    normalized[face.getIndex()] = Vec3(faceCoords.x, faceCoords.y, faceCoords.z);
  }
  geometry.unrequireEdgeIndices();
  return normalized;
}

std::vector<Vec3> averageFaceDirectionsToVerticesReference(
    GeometryCentralSurface& surface,
    const FaceHeatDirectionField& normalizedFaceDirections) {
  auto& mesh = *surface.mesh;
  auto& geometry = *surface.geometry;
  if (normalizedFaceDirections.size() != mesh.nFaces()) {
    throw std::runtime_error("averageFaceDirectionsToVerticesReference requires one direction per face");
  }

  geometry.requireFaceAreas();
  geometry.requireVertexPositions();
  std::vector<Vec3> vertexVectors(mesh.nVertices(), Vec3::Zero());
  for (gcs::Vertex vertex : mesh.vertices()) {
    gc::Vector3 average{0.0, 0.0, 0.0};
    double totalArea = 0.0;
    for (gcs::Face face : vertex.adjacentFaces()) {
      const Vec3& coords = normalizedFaceDirections[face.getIndex()];
      const gcs::SurfacePoint point(face, gc::Vector3{coords.x(), coords.y(), coords.z()});
      const double area = geometry.faceAreas[face];
      average += area * point.interpolate(geometry.vertexPositions);
      totalArea += area;
    }
    if (totalArea > 0.0) {
      average /= totalArea;
    }
    vertexVectors[vertex.getIndex()] = fromGcVector3(average);
  }
  geometry.unrequireFaceAreas();
  geometry.unrequireVertexPositions();
  return vertexVectors;
}

std::vector<Vec3> averageFaceDirectionsToVerticesProjected(
    GeometryCentralSurface& surface,
    const FaceHeatDirectionField& normalizedFaceDirections,
    VertexDirectionAveraging averaging) {
  auto& mesh = *surface.mesh;
  auto& geometry = *surface.geometry;
  if (normalizedFaceDirections.size() != mesh.nFaces()) {
    throw std::runtime_error("averageFaceDirectionsToVerticesProjected requires one direction per face");
  }

  geometry.requireFaceAreas();
  geometry.requireVertexPositions();

  std::vector<gc::Vector3> accumulated(mesh.nVertices(), gc::Vector3{0.0, 0.0, 0.0});
  std::vector<gc::Vector3> accumulatedNormals(mesh.nVertices(), gc::Vector3{0.0, 0.0, 0.0});
  std::vector<double> totalWeights(mesh.nVertices(), 0.0);

  for (gcs::Face face : mesh.faces()) {
    const Vec3& coords = normalizedFaceDirections[face.getIndex()];
    const gc::Vector3 faceVector = faceVectorToExtrinsic(face, coords, geometry);

    std::array<gcs::Vertex, 3> vertices;
    std::array<gc::Vector3, 3> positions;
    size_t localIndex = 0;
    for (gcs::Vertex vertex : face.adjacentVertices()) {
      vertices[localIndex] = vertex;
      positions[localIndex] = geometry.vertexPositions[vertex];
      ++localIndex;
    }

    const gc::Vector3 faceNormal = cross(positions[1] - positions[0], positions[2] - positions[0]);
    const double faceArea = 0.5 * norm(faceNormal);
    if (faceArea == 0.0) {
      continue;
    }

    for (size_t i = 0; i < vertices.size(); ++i) {
      double weight = faceArea;
      if (averaging == VertexDirectionAveraging::CornerAngle) {
        weight = cornerAngle(positions[i], positions[(i + 1) % 3], positions[(i + 2) % 3]);
      }
      const size_t vertexIndex = vertices[i].getIndex();
      accumulated[vertexIndex] += weight * faceVector;
      accumulatedNormals[vertexIndex] += faceNormal;
      totalWeights[vertexIndex] += weight;
    }
  }

  std::vector<Vec3> vertexVectors(mesh.nVertices(), Vec3::Zero());
  for (gcs::Vertex vertex : mesh.vertices()) {
    const size_t vertexIndex = vertex.getIndex();
    if (totalWeights[vertexIndex] == 0.0) {
      continue;
    }

    gc::Vector3 vector = accumulated[vertexIndex] / totalWeights[vertexIndex];
    const double normalNorm = norm(accumulatedNormals[vertexIndex]);
    if (normalNorm > 0.0) {
      const gc::Vector3 normal = accumulatedNormals[vertexIndex] / normalNorm;
      vector -= dot(vector, normal) * normal;
    }
    const double vectorNorm = norm(vector);
    if (vectorNorm > 0.0) {
      vector /= vectorNorm;
    }
    vertexVectors[vertexIndex] = fromGcVector3(vector);
  }

  geometry.unrequireFaceAreas();
  geometry.unrequireVertexPositions();
  return vertexVectors;
}

CustomSignedHeatResult computeCustomSignedHeatDirections(GeometryCentralSurface& surface,
                                               const std::vector<SurfaceReference>& sourceCurve,
                                               const SignedHeatSolveOptions& options) {
  CustomSignedHeatSolver solver(surface, options.diffusionTimeCoefficient);
  CustomSignedHeatResult result;
  result.diffusion = solver.solveDiffusedEdgeHeatField(sourceCurve, options);
  result.normalizedFaceDirections =
      sampleAndNormalizeFaceDirections(surface, result.diffusion.diffusedEdgeHeatField);
  result.vertexDirections =
      averageFaceDirectionsToVerticesReference(surface, result.normalizedFaceDirections);
  return result;
}

std::array<CustomSignedHeatResult, 2> computeCustomSignedHeatDirections(GeometryCentralSurface& surface,
                                                                        const SourceCurves& sourceCurves,
                                                                        const SignedHeatSolveOptions& options) {
  CustomSignedHeatSolver solver(surface, options.diffusionTimeCoefficient);
  std::array<CustomSignedHeatResult, 2> results;
  for (size_t i = 0; i < results.size(); ++i) {
    results[i].diffusion = solver.solveDiffusedEdgeHeatField(sourceCurves.curves[i], options);
    results[i].normalizedFaceDirections =
        sampleAndNormalizeFaceDirections(surface, results[i].diffusion.diffusedEdgeHeatField);
    results[i].vertexDirections =
        averageFaceDirectionsToVerticesReference(surface, results[i].normalizedFaceDirections);
  }
  return results;
}

} // namespace geodesic_draping


#include "geodrape_internal.h"

#include "geometrycentral/surface/integer_coordinates_intrinsic_triangulation.h"
#include "geometrycentral/surface/signpost_intrinsic_triangulation.h"

#include <cmath>
#include <limits>
#include <stdexcept>

namespace geodesic_draping {
namespace {

std::unique_ptr<gcs::IntrinsicTriangulation> makeIntrinsicTriangulation(
    IntrinsicTriangulationBackend backend,
    gcs::ManifoldSurfaceMesh& mesh,
    gcs::IntrinsicGeometryInterface& geometry) {
  if (backend == IntrinsicTriangulationBackend::Signpost) {
    return std::make_unique<gcs::SignpostIntrinsicTriangulation>(mesh, geometry);
  }
  return std::make_unique<gcs::IntegerCoordinatesIntrinsicTriangulation>(mesh, geometry);
}

Vec3 toEigenVector(const geometrycentral::Vector3& value) {
  return Vec3(value.x, value.y, value.z);
}

std::optional<Vec3> barycentricInBarycentricTriangle(const Vec3& point,
                                                     const Vec3& a,
                                                     const Vec3& b,
                                                     const Vec3& c) {
  Eigen::Matrix3d matrix;
  matrix.col(0) = a;
  matrix.col(1) = b;
  matrix.col(2) = c;
  const double determinant = matrix.determinant();
  if (std::abs(determinant) < 1e-14) {
    return std::nullopt;
  }
  return matrix.inverse() * point;
}

gcs::SurfacePoint inputToIntrinsicViaCommonSubdivision(gcs::IntrinsicTriangulation& triangulation,
                                                       const gcs::SurfacePoint& pointOnInput) {
  const gcs::SurfacePoint inputFacePoint = pointOnInput.inSomeFace();
  gcs::CommonSubdivision& subdivision = triangulation.getCommonSubdivision();
  subdivision.constructMesh();

  const Vec3 query = toEigenVector(inputFacePoint.faceCoords);
  for (gcs::Face subdivisionFace : subdivision.mesh->faces()) {
    if (subdivision.sourceFaceA[subdivisionFace] != inputFacePoint.face) {
      continue;
    }

    const gcs::Face intrinsicFace = subdivision.sourceFaceB[subdivisionFace];
    std::array<Vec3, 3> pointsA{};
    std::array<Vec3, 3> pointsB{};
    size_t i = 0;
    for (gcs::Vertex vertex : subdivisionFace.adjacentVertices()) {
      const gcs::CommonSubdivisionPoint* source = subdivision.sourcePoints[vertex];
      pointsA[i] = toEigenVector(source->posA.inFace(inputFacePoint.face).faceCoords);
      pointsB[i] = toEigenVector(source->posB.inFace(intrinsicFace).faceCoords);
      ++i;
    }

    const std::optional<Vec3> localBary =
        barycentricInBarycentricTriangle(query, pointsA[0], pointsA[1], pointsA[2]);
    if (!localBary) {
      continue;
    }

    constexpr double tolerance = 1e-8;
    if ((*localBary).minCoeff() < -tolerance || (*localBary).maxCoeff() > 1.0 + tolerance) {
      continue;
    }

    Vec3 intrinsicBary =
        (*localBary).x() * pointsB[0] +
        (*localBary).y() * pointsB[1] +
        (*localBary).z() * pointsB[2];
    intrinsicBary = intrinsicBary.cwiseMax(0.0);
    const double sum = intrinsicBary.sum();
    if (sum <= 0.0) {
      continue;
    }
    intrinsicBary /= sum;
    return gcs::SurfacePoint(
        intrinsicFace,
        geometrycentral::Vector3{intrinsicBary.x(), intrinsicBary.y(), intrinsicBary.z()});
  }

  throw std::runtime_error("failed to map input point to intrinsic triangulation via common subdivision");
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
                                             const IntrinsicConstructionOptions& intrinsicOptions,
                                             const RefinementOptions& refinementOptions) {
  useCommonSubdivisionInputAdapter_ =
      intrinsicOptions.backend == IntrinsicTriangulationBackend::IntegerCoordinates;
  triangulation_ = makeIntrinsicTriangulation(
      intrinsicOptions.backend,
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
  if (useCommonSubdivisionInputAdapter_) {
    return inputToIntrinsicViaCommonSubdivision(*triangulation_, pointOnInput);
  }
  return triangulation_->equivalentPointOnIntrinsic(pointOnInput);
}

gcs::SurfacePoint ActiveIntrinsicDomain::intrinsicToInput(const gcs::SurfacePoint& pointOnIntrinsic) {
  return triangulation_->equivalentPointOnInput(pointOnIntrinsic);
}

} // namespace geodesic_draping

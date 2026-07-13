#include "geodrape_internal.h"

#include "geometrycentral/surface/integer_coordinates_intrinsic_triangulation.h"
#include "geometrycentral/surface/signpost_intrinsic_triangulation.h"

#include <cmath>
#include <limits>
#include <stdexcept>

namespace geodesic_draping {
namespace {

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

std::unique_ptr<gcs::IntrinsicTriangulation> makeIntrinsicTriangulation(
    IntrinsicTriangulationBackend backend,
    gcs::ManifoldSurfaceMesh& mesh,
    gcs::IntrinsicGeometryInterface& geometry) {
  if (backend == IntrinsicTriangulationBackend::Signpost) {
    return std::make_unique<gcs::SignpostIntrinsicTriangulation>(mesh, geometry);
  }
  return std::make_unique<gcs::IntegerCoordinatesIntrinsicTriangulation>(mesh, geometry);
}

void applyRefinement(gcs::IntrinsicTriangulation& triangulation,
                     const RefinementOptions& refinementOptions) {
  if (refinementOptions.mode == RefinementMode::DelaunayFlip) {
    triangulation.flipToDelaunay();
  } else if (refinementOptions.mode == RefinementMode::DelaunayRefine) {
    const double angleThreshold = refinementOptions.angleThreshold.value_or(25.0);
    const double circumradiusThreshold =
        refinementOptions.circumradiusThreshold.value_or(std::numeric_limits<double>::infinity());
    const size_t maxInsertions =
        refinementOptions.maxInsertions.value_or(geometrycentral::INVALID_IND);
    triangulation.delaunayRefine(angleThreshold, circumradiusThreshold, maxInsertions);
  }
}

gcs::SurfacePoint inputToIntrinsic(gcs::IntrinsicTriangulation& triangulation,
                                   const gcs::SurfacePoint& pointOnInput,
                                   bool useCommonSubdivisionInputAdapter) {
  if (useCommonSubdivisionInputAdapter) {
    // geometry-central's IntegerCoordinatesIntrinsicTriangulation currently does not
    // support equivalentPointOnIntrinsic(), so input points are mapped through the
    // common subdivision instead.
    return inputToIntrinsicViaCommonSubdivision(triangulation, pointOnInput);
  }
  return triangulation.equivalentPointOnIntrinsic(pointOnInput);
}

} // namespace geodesic_draping

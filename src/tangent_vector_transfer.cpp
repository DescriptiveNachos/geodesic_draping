#include "geodrape_internal.h"

#include <stdexcept>

namespace geodesic_draping {

gcs::BarycentricVector toBarycentricVector(gcs::SurfaceMesh& mesh,
                                           const TangentVectorRef& ref) {
  if (ref.type != SurfaceReferenceType::Face || ref.coords.size() != 3) {
    throw std::runtime_error("solveFromIntrinsic currently requires a face tangent vector with three barycentric coordinates");
  }
  return gcs::BarycentricVector(
      mesh.face(ref.elementIndex),
      geometrycentral::Vector3{ref.coords[0], ref.coords[1], ref.coords[2]});
}

TangentVectorRef toTangentVectorRef(const gcs::BarycentricVector& vector) {
  const gcs::BarycentricVector faceVector = vector.inSomeFace();
  TangentVectorRef ref;
  ref.type = SurfaceReferenceType::Face;
  ref.elementIndex = faceVector.face.getIndex();
  ref.coords = {faceVector.faceCoords.x, faceVector.faceCoords.y, faceVector.faceCoords.z};
  return ref;
}

std::array<Vec3, 4> cartesianDirectionsFromIntrinsic(
    gcs::VertexPositionGeometry& inputGeometry,
    gcs::IntrinsicTriangulation& triangulation,
    const gcs::SurfacePoint& seed,
    const std::array<gcs::BarycentricVector, 4>& directions) {
  std::array<Vec3, 4> out{};
  const gcs::SurfacePoint seedFace = seed.inSomeFace();
  std::array<Vec3, 3> positions;
  size_t vertexIndex = 0;
  for (gcs::Vertex vertex : seedFace.face.adjacentVertices()) {
    const gcs::SurfacePoint inputPoint = triangulation.equivalentPointOnInput(gcs::SurfacePoint(vertex));
    positions[vertexIndex] = interpolateSurfacePoint(inputPoint, inputGeometry);
    ++vertexIndex;
  }

  for (size_t i = 0; i < directions.size(); ++i) {
    const gcs::BarycentricVector inFace = directions[i].inFace(seedFace.face);
    Vec3 delta = inFace.faceCoords.x * positions[0] +
                 inFace.faceCoords.y * positions[1] +
                 inFace.faceCoords.z * positions[2];
    const double norm = delta.norm();
    out[i] = Vec3::Zero();
    if (norm > 0.0) {
      out[i] = delta / norm;
    }
  }
  return out;
}

} // namespace geodesic_draping

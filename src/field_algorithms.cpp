#include "geodrape_internal.h"

#include <stdexcept>

namespace geodesic_draping {

FaceHeatDirectionField computeIntrinsicFaceScalarGradients(gcs::SurfaceMesh& mesh,
                                                           gcs::IntrinsicGeometryInterface& geometry,
                                                           const std::vector<double>& scalarField) {
  if (scalarField.size() != mesh.nVertices()) {
    throw std::runtime_error("computeIntrinsicFaceScalarGradients requires one scalar per active vertex");
  }

  FaceHeatDirectionField gradients(mesh.nFaces(), geometrycentral::Vector3{0.0, 0.0, 0.0});
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
    gradients[face.getIndex()] = gradient.faceCoords;
  }
  return gradients;
}

} // namespace geodesic_draping

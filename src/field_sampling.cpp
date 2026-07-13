#include "geodrape_internal.h"

#include <stdexcept>

namespace geodesic_draping {

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

} // namespace geodesic_draping

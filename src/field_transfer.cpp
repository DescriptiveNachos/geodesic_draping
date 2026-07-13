#include "geodrape_internal.h"

#include <stdexcept>

namespace geodesic_draping {

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

FaceHeatDirectionField faceVectorDataToVector(const gcs::FaceData<geometrycentral::Vector3>& values) {
  FaceHeatDirectionField out(values.getMesh()->nFaces(), geometrycentral::Vector3{0.0, 0.0, 0.0});
  for (gcs::Face face : values.getMesh()->faces()) {
    out[face.getIndex()] = values[face];
  }
  return out;
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
  return vertexDataToVector(triangulation.restrictToInput(activeValues));
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

gcs::FaceData<geometrycentral::Vector3> activeFaceVectorData(gcs::SurfaceMesh& mesh,
                                                             const FaceHeatDirectionField& values) {
  if (values.size() != mesh.nFaces()) {
    throw std::runtime_error("activeFaceVectorData requires one vector per active face");
  }
  gcs::FaceData<geometrycentral::Vector3> data(mesh, geometrycentral::Vector3{0.0, 0.0, 0.0});
  for (gcs::Face face : mesh.faces()) {
    data[face] = values[face.getIndex()];
  }
  return data;
}

} // namespace geodesic_draping

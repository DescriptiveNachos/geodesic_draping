#include "geodrape_internal.h"

#include <stdexcept>

namespace geodesic_draping {

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

} // namespace geodesic_draping

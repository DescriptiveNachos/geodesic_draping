#include "geodesic_draping/field_processing.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace geodesic_draping {

std::vector<Vec3> computeVertexScalarGradients(const SurfaceMeshData& mesh,
                                               const std::vector<double>& scalarField) {
  if (mesh.vertices.size() != scalarField.size()) {
    throw std::runtime_error("computeVertexScalarGradients requires one scalar per vertex");
  }

  std::vector<Vec3> accumulated(mesh.vertices.size(), Vec3::Zero());
  std::vector<double> areaWeights(mesh.vertices.size(), 0.0);

  for (const Face& face : mesh.faces) {
    const Vec3& p0 = mesh.vertices[face[0]];
    const Vec3& p1 = mesh.vertices[face[1]];
    const Vec3& p2 = mesh.vertices[face[2]];
    const Vec3 e1 = p1 - p0;
    const Vec3 e2 = p2 - p0;
    const Vec3 normalCross = e1.cross(e2);
    const double area = 0.5 * normalCross.norm();
    if (area == 0.0) {
      continue;
    }

    Eigen::Matrix2d gram;
    gram << e1.dot(e1), e1.dot(e2),
            e2.dot(e1), e2.dot(e2);
    Eigen::Vector2d rhs;
    rhs << scalarField[face[1]] - scalarField[face[0]],
           scalarField[face[2]] - scalarField[face[0]];
    const Eigen::Vector2d uv = gram.ldlt().solve(rhs);
    const Vec3 faceGradient = uv.x() * e1 + uv.y() * e2;

    for (size_t vertexIndex : face) {
      accumulated[vertexIndex] += area * faceGradient;
      areaWeights[vertexIndex] += area;
    }
  }

  for (size_t i = 0; i < accumulated.size(); ++i) {
    if (areaWeights[i] > 0.0) {
      accumulated[i] /= areaWeights[i];
    }
  }
  return accumulated;
}

std::vector<double> computeShearAnglesDegrees(const std::vector<Vec3>& gradients0,
                                              const std::vector<Vec3>& gradients1) {
  if (gradients0.size() != gradients1.size()) {
    throw std::runtime_error("computeShearAnglesDegrees requires equally sized gradient arrays");
  }

  constexpr double pi = 3.141592653589793238462643383279502884;
  std::vector<double> shear;
  shear.reserve(gradients0.size());
  for (size_t i = 0; i < gradients0.size(); ++i) {
    const double norm0 = gradients0[i].norm();
    const double norm1 = gradients1[i].norm();
    const double cosTheta = std::clamp(gradients1[i].dot(gradients0[i]) / (norm0 * norm1), -1.0, 1.0);
    shear.push_back(std::abs((std::acos(cosTheta) * 180.0 / pi) - 90.0));
  }
  return shear;
}

} // namespace geodesic_draping

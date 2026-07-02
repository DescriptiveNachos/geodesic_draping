#include "geodesic_draping/field_processing.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace geodesic_draping {
namespace {

double triangleArea(const Vec3& p0, const Vec3& p1, const Vec3& p2) {
  return 0.5 * (p1 - p0).cross(p2 - p0).norm();
}

double cornerAngle(const Vec3& center, const Vec3& pA, const Vec3& pB) {
  const Vec3 u = pA - center;
  const Vec3 v = pB - center;
  const double denom = u.norm() * v.norm();
  if (denom == 0.0) {
    return 0.0;
  }
  return std::acos(std::clamp(u.dot(v) / denom, -1.0, 1.0));
}

} // namespace

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
    if (norm0 == 0.0 || norm1 == 0.0) {
      shear.push_back(std::numeric_limits<double>::quiet_NaN());
      continue;
    }
    const double dotProduct = gradients0[i].dot(gradients1[i]);
    const double determinantMagnitude = gradients0[i].cross(gradients1[i]).norm();
    const double theta = std::atan2(determinantMagnitude, dotProduct);
    shear.push_back(std::abs((theta * 180.0 / pi) - 90.0));
  }
  return shear;
}

std::vector<double> averageFaceScalarsToVertices(const SurfaceMeshData& mesh,
                                                 const std::vector<double>& faceScalars,
                                                 FaceScalarAveraging averaging) {
  if (mesh.faces.size() != faceScalars.size()) {
    throw std::runtime_error("averageFaceScalarsToVertices requires one scalar per face");
  }

  std::vector<double> accumulated(mesh.vertices.size(), 0.0);
  std::vector<double> totalWeights(mesh.vertices.size(), 0.0);
  for (size_t faceIndex = 0; faceIndex < mesh.faces.size(); ++faceIndex) {
    const Face& face = mesh.faces[faceIndex];
    const Vec3& p0 = mesh.vertices[face[0]];
    const Vec3& p1 = mesh.vertices[face[1]];
    const Vec3& p2 = mesh.vertices[face[2]];
    const double area = triangleArea(p0, p1, p2);
    if (area == 0.0) {
      continue;
    }

    const std::array<double, 3> cornerWeights{
        averaging == FaceScalarAveraging::CornerAngle ? cornerAngle(p0, p1, p2) : area,
        averaging == FaceScalarAveraging::CornerAngle ? cornerAngle(p1, p2, p0) : area,
        averaging == FaceScalarAveraging::CornerAngle ? cornerAngle(p2, p0, p1) : area,
    };
    for (size_t localIndex = 0; localIndex < 3; ++localIndex) {
      const size_t vertexIndex = face[localIndex];
      const double weight = cornerWeights[localIndex];
      accumulated[vertexIndex] += weight * faceScalars[faceIndex];
      totalWeights[vertexIndex] += weight;
    }
  }

  for (size_t i = 0; i < accumulated.size(); ++i) {
    if (totalWeights[i] > 0.0) {
      accumulated[i] /= totalWeights[i];
    }
  }
  return accumulated;
}

std::vector<double> averageVertexScalarsToFaces(const SurfaceMeshData& mesh,
                                                const std::vector<double>& vertexScalars) {
  if (mesh.vertices.size() != vertexScalars.size()) {
    throw std::runtime_error("averageVertexScalarsToFaces requires one scalar per vertex");
  }

  std::vector<double> faceScalars;
  faceScalars.reserve(mesh.faces.size());
  for (const Face& face : mesh.faces) {
    faceScalars.push_back(
        (vertexScalars[face[0]] + vertexScalars[face[1]] + vertexScalars[face[2]]) / 3.0);
  }
  return faceScalars;
}

} // namespace geodesic_draping

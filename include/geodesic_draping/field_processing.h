#pragma once

#include "geodesic_draping/mesh.h"

#include <vector>

namespace geodesic_draping {

enum class FaceScalarAveraging {
  FaceArea,
  CornerAngle,
};

std::vector<Vec3> computeVertexScalarGradients(const SurfaceMeshData& mesh,
                                               const std::vector<double>& scalarField);

std::vector<double> computeShearAnglesDegrees(const std::vector<Vec3>& gradients0,
                                              const std::vector<Vec3>& gradients1);

std::vector<double> averageFaceScalarsToVertices(const SurfaceMeshData& mesh,
                                                 const std::vector<double>& faceScalars,
                                                 FaceScalarAveraging averaging = FaceScalarAveraging::FaceArea);

// Samples a vertex scalar field at face centroids. This is a smoothing projection,
// not an inverse of averageFaceScalarsToVertices().
std::vector<double> averageVertexScalarsToFaces(const SurfaceMeshData& mesh,
                                                const std::vector<double>& vertexScalars);

} // namespace geodesic_draping

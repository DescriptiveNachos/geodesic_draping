#pragma once

#include "geodesic_draping/geodrape.h"

#include "geometrycentral/surface/common_subdivision.h"

#include <array>
#include <memory>
#include <optional>
#include <vector>

namespace geodesic_draping {
namespace gcs = geometrycentral::surface;

inline constexpr double kPi = 3.141592653589793238462643383279502884;

TraceSettings resolveTraceSettings(const TraceSettings& defaults,
                                   const AdvancedTraceOptions& overrides);
double boundingBoxDiagonal(const SurfaceMeshData& meshData);
TraceSettings makeTraceDefaults(const SurfaceMeshData& meshData,
                                const gcs::SurfaceMesh& activeMesh);

std::unique_ptr<gcs::IntrinsicTriangulation> makeIntrinsicTriangulation(
    IntrinsicTriangulationBackend backend,
    gcs::ManifoldSurfaceMesh& mesh,
    gcs::IntrinsicGeometryInterface& geometry);
void applyRefinement(gcs::IntrinsicTriangulation& triangulation,
                     const RefinementOptions& refinementOptions);
gcs::SurfacePoint inputToIntrinsic(
    gcs::IntrinsicTriangulation& triangulation,
    const gcs::SurfacePoint& pointOnInput,
    bool useCommonSubdivisionInputAdapter);

gcs::SurfacePoint toFaceSurfacePoint(gcs::SurfaceMesh& mesh, const BarycentricPoint& point);
Vec3 interpolateSurfacePoint(const gcs::SurfacePoint& point,
                             gcs::VertexPositionGeometry& geometry);

gcs::BarycentricVector normalizeVector(gcs::BarycentricVector vector,
                                       gcs::IntrinsicGeometryInterface& geometry);
gcs::BarycentricVector intrinsicDirectionFromFabricAngle(const SurfaceMeshData& meshData,
                                                        gcs::VertexPositionGeometry& inputGeometry,
                                                        gcs::IntrinsicTriangulation& triangulation,
                                                        const gcs::SurfacePoint& inputSeed,
                                                        const gcs::SurfacePoint& intrinsicSeed,
                                                        double fabricAngle,
                                                        bool inputConnectivityPreserved,
                                                        bool useCommonSubdivisionInputAdapter);

std::array<IntrinsicGeneratorTrace, 4> traceIntrinsicGenerators(
    gcs::IntrinsicTriangulation& triangulation,
    const gcs::SurfacePoint& start,
    const std::array<gcs::BarycentricVector, 4>& directions,
    const TraceSettings& settings);
std::array<gcs::Curve, 2> pairOppositeIntrinsicGeneratorTraces(
    const std::array<IntrinsicGeneratorTrace, 4>& traces);

// Core field algorithms.
FaceHeatDirectionField computeIntrinsicFaceScalarGradients(
    gcs::SurfaceMesh& mesh,
    gcs::IntrinsicGeometryInterface& geometry,
    const std::vector<double>& scalarField);

// Retrieval-time sampling.
std::vector<double> averageIntrinsicFaceScalarsToVertices(
    gcs::SurfaceMesh& mesh,
    const std::vector<double>& faceScalars);

} // namespace geodesic_draping

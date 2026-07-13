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

gcs::BarycentricVector normalizeVector(gcs::BarycentricVector vector,
                                       gcs::IntrinsicGeometryInterface& geometry);
gcs::BarycentricVector toBarycentricVector(gcs::SurfaceMesh& mesh,
                                           const TangentVectorRef& ref);
TangentVectorRef toTangentVectorRef(const gcs::BarycentricVector& vector);
std::array<Vec3, 4> cartesianDirectionsFromIntrinsic(
    gcs::VertexPositionGeometry& inputGeometry,
    gcs::IntrinsicTriangulation& triangulation,
    const gcs::SurfacePoint& seed,
    const std::array<gcs::BarycentricVector, 4>& directions);
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

// API/result transfer helpers.
SeedProjection toExtrinsicSeed(gcs::VertexPositionGeometry& inputGeometry,
                               gcs::IntrinsicTriangulation& triangulation,
                               const gcs::SurfacePoint& intrinsicSeed);
std::vector<double> restrictVertexScalarsToInput(
    gcs::IntrinsicTriangulation& triangulation,
    const std::vector<double>& valuesOnIntrinsic);

ResultMesh makeExtrinsicResultMesh(const SurfaceMeshData& meshData);
ResultMesh makeIntrinsicResultMesh(gcs::SurfaceMesh& mesh,
                                   gcs::IntrinsicGeometryInterface& geometry);
ResultMesh makeSubdivisionResultMesh(gcs::CommonSubdivision& subdivision,
                                     gcs::VertexPositionGeometry& inputGeometry);
std::array<TraceFamily, 2> makeTraceFamilies(
    gcs::VertexPositionGeometry& inputGeometry,
    gcs::IntrinsicTriangulation& triangulation,
    const std::array<IntrinsicGeneratorTrace, 4>& generators,
    ResultDomain domain);

std::vector<double> vertexDataToVector(const gcs::VertexData<double>& values);
std::vector<double> faceDataToVector(const gcs::FaceData<double>& values);
FaceHeatDirectionField faceVectorDataToVector(const gcs::FaceData<geometrycentral::Vector3>& values);
gcs::VertexData<double> activeVertexData(gcs::SurfaceMesh& mesh, const std::vector<double>& values);
gcs::FaceData<double> activeFaceData(gcs::SurfaceMesh& mesh, const std::vector<double>& values);
gcs::FaceData<geometrycentral::Vector3> activeFaceVectorData(gcs::SurfaceMesh& mesh,
                                                             const FaceHeatDirectionField& values);

} // namespace geodesic_draping

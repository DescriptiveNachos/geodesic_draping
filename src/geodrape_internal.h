#pragma once

#include "geodesic_draping/geodrape.h"

#include "geometrycentral/surface/common_subdivision.h"

#include <array>
#include <optional>
#include <vector>

namespace geodesic_draping {
namespace gcs = geometrycentral::surface;

inline constexpr double kPi = 3.141592653589793238462643383279502884;

std::array<Vec3, 4> generateCartesianFamilyDirections(double fabricAngle,
                                                       double fiberAngle);
TraceSettings resolveTraceSettings(const TraceSettings& defaults,
                                   const AdvancedTraceOptions& overrides);
double boundingBoxDiagonal(const SurfaceMeshData& meshData);
TraceSettings makeTraceDefaults(const SurfaceMeshData& meshData,
                                const gcs::SurfaceMesh& activeMesh);

gcs::SurfacePoint toFaceSurfacePoint(gcs::SurfaceMesh& mesh, const BarycentricPoint& point);
SeedProjection toExtrinsicSeed(ReferenceGeometry& reference,
                               ActiveIntrinsicDomain& activeDomain,
                               const SurfaceReference& intrinsicSeed);

gcs::BarycentricVector normalizeVector(gcs::BarycentricVector vector,
                                       gcs::IntrinsicGeometryInterface& geometry);
gcs::BarycentricVector toBarycentricVector(gcs::SurfaceMesh& mesh,
                                           const TangentVectorRef& ref);
TangentVectorRef toTangentVectorRef(const gcs::BarycentricVector& vector);
std::array<Vec3, 4> cartesianDirectionsFromIntrinsic(
    ReferenceGeometry& reference,
    ActiveIntrinsicDomain& activeDomain,
    const gcs::SurfacePoint& seed,
    const std::array<gcs::BarycentricVector, 4>& directions);
gcs::BarycentricVector intrinsicDirectionFromFabricAngle(ReferenceGeometry& reference,
                                                        ActiveIntrinsicDomain& activeDomain,
                                                        const gcs::SurfacePoint& inputSeed,
                                                        const gcs::SurfacePoint& intrinsicSeed,
                                                        double fabricAngle,
                                                        bool inputConnectivityPreserved);

std::array<GeneratorTrace, 4> traceActiveGenerators(
    ReferenceGeometry& reference,
    ActiveIntrinsicDomain& activeDomain,
    const gcs::SurfacePoint& start,
    const std::array<gcs::BarycentricVector, 4>& directions,
    const TraceSettings& settings);

FaceHeatDirectionField computeIntrinsicFaceScalarGradients(
    gcs::SurfaceMesh& mesh,
    gcs::IntrinsicGeometryInterface& geometry,
    const std::vector<double>& scalarField);
std::vector<double> averageIntrinsicFaceScalarsToVertices(
    gcs::SurfaceMesh& mesh,
    const std::vector<double>& faceScalars);
std::vector<double> restrictVertexScalarsToInput(
    gcs::IntrinsicTriangulation& triangulation,
    const std::vector<double>& valuesOnIntrinsic);

ResultMesh makeExtrinsicResultMesh(const SurfaceMeshData& meshData);
ResultMesh makeIntrinsicResultMesh(gcs::SurfaceMesh& mesh,
                                   gcs::IntrinsicGeometryInterface& geometry);
ResultMesh makeSubdivisionResultMesh(gcs::CommonSubdivision& subdivision,
                                     gcs::VertexPositionGeometry& inputGeometry);
std::array<TraceFamily, 2> makeTraceFamilies(const std::array<GeneratorTrace, 4>& generators,
                                             ResultDomain domain);

std::vector<double> vertexDataToVector(const gcs::VertexData<double>& values);
std::vector<double> faceDataToVector(const gcs::FaceData<double>& values);
FaceHeatDirectionField faceVectorDataToVector(const gcs::FaceData<Vec3>& values);
gcs::VertexData<double> activeVertexData(gcs::SurfaceMesh& mesh, const std::vector<double>& values);
gcs::FaceData<double> activeFaceData(gcs::SurfaceMesh& mesh, const std::vector<double>& values);
gcs::FaceData<Vec3> activeFaceVectorData(gcs::SurfaceMesh& mesh, const FaceHeatDirectionField& values);

} // namespace geodesic_draping

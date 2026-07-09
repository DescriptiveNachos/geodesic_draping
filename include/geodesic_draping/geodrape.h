#pragma once

#include "geodesic_draping/field_processing.h"
#include "geodesic_draping/generator_tracing.h"
#include "geodesic_draping/custom_signed_heat.h"

#include "geometrycentral/surface/intrinsic_triangulation.h"
#include "geometrycentral/surface/barycentric_vector.h"
#include "geometrycentral/surface/surface_point.h"

#include <array>
#include <memory>
#include <optional>
#include <string>
#include <vector>

namespace geodesic_draping {

enum class DrapeSolveMode {
  Fast,
  Hybrid,
  Complete,
};

enum class ResultDomain {
  Intrinsic,
  Extrinsic,
  Subdivision,
};

enum class RefinementMode {
  None,
  DelaunayFlip,
  DelaunayRefine,
};

struct IntrinsicConstructionOptions {};

struct RefinementOptions {
  RefinementMode mode = RefinementMode::None;
  std::optional<double> angleThreshold;
  std::optional<double> circumradiusThreshold;
  std::optional<size_t> maxInsertions;
};

struct AdvancedTraceOptions {
  std::optional<double> traceLength;
  std::optional<size_t> maxIterations;
};

struct AdvancedSolveOptions {
  AdvancedTraceOptions trace;
};

struct DrapeSolveOptions {
  DrapeSolveMode mode = DrapeSolveMode::Complete;
  ResultDomain retrieval = ResultDomain::Extrinsic;
  bool sampleVertexShear = true;
  AdvancedSolveOptions advanced;
};

using GluingMapEntry = std::array<int, 2>;
using FaceGluingMap = std::array<GluingMapEntry, 3>;

struct ResultMesh {
  ResultDomain domain = ResultDomain::Extrinsic;
  std::vector<Face> faces;
  std::optional<std::vector<Vec3>> vertices3D;
  std::optional<std::vector<std::array<double, 3>>> edgeLengths;
  std::optional<std::vector<FaceGluingMap>> gluingMap;
  std::optional<std::vector<double>> intrinsicSourceFaceColor;
};

struct TangentVectorRef {
  SurfaceReferenceType type = SurfaceReferenceType::Face;
  size_t elementIndex = 0;
  std::vector<double> coords;
};

struct DrapeOrigin {
  std::optional<SurfaceReference> intrinsicPoint;
  std::optional<Vec3> extrinsicPoint;
  std::array<std::optional<TangentVectorRef>, 2> intrinsicFamilyDirections;
  std::array<std::optional<Vec3>, 2> extrinsicFamilyDirections;
};

struct DrapeTrace {
  std::vector<SurfaceReference> intrinsicPoints;
  std::vector<Vec3> extrinsicPoints;
  bool hitBoundary = false;
  double length = 0.0;
};

struct TraceFamily {
  DrapeTrace positive;
  DrapeTrace negative;
};

struct DrapeResult {
  ResultDomain domain = ResultDomain::Extrinsic;
  DrapeSolveMode mode = DrapeSolveMode::Fast;
  ResultMesh mesh;
  DrapeOrigin origin;
  std::array<TraceFamily, 2> traces;

  std::array<FaceHeatDirectionField, 2> directions;
  std::optional<std::array<std::vector<double>, 2>> distances;
  std::optional<std::vector<double>> faceShear;
  std::optional<std::vector<double>> vertexShear;
};

class ReferenceGeometry {
public:
  explicit ReferenceGeometry(SurfaceMeshData meshData);

  SurfaceMeshData& meshData();
  const SurfaceMeshData& meshData() const;
  GeometryCentralSurface& surface();
  const GeometryCentralSurface& surface() const;

private:
  SurfaceMeshData meshData_;
  GeometryCentralSurface surface_;
};

class ActiveIntrinsicDomain {
public:
  ActiveIntrinsicDomain(ReferenceGeometry& reference,
                        const IntrinsicConstructionOptions& intrinsicOptions,
                        const RefinementOptions& refinementOptions);

  geometrycentral::surface::ManifoldSurfaceMesh& mesh();
  geometrycentral::surface::IntrinsicGeometryInterface& geometry();
  geometrycentral::surface::IntrinsicTriangulation& triangulation();
  geometrycentral::surface::SurfacePoint inputToIntrinsic(
      const geometrycentral::surface::SurfacePoint& pointOnInput);
  geometrycentral::surface::SurfacePoint intrinsicToInput(
      const geometrycentral::surface::SurfacePoint& pointOnIntrinsic);

private:
  std::unique_ptr<geometrycentral::surface::IntrinsicTriangulation> triangulation_;
};

struct IntrinsicSolveInput {
  geometrycentral::surface::SurfacePoint seed;
  std::array<geometrycentral::surface::BarycentricVector, 4> directions;
  std::array<Vec3, 4> cartesianDirections;
  DrapeSolveMode mode = DrapeSolveMode::Complete;
  TraceSettings trace;
};

struct CoreIntrinsicResult {
  DrapeSolveMode mode = DrapeSolveMode::Complete;
  SurfaceReference intrinsicSeed;
  std::array<TangentVectorRef, 4> intrinsicDirections;
  std::array<Vec3, 4> cartesianDirections;
  std::array<GeneratorTrace, 4> generators;
  std::array<FaceHeatDirectionField, 2> directions;
  std::optional<std::array<std::vector<double>, 2>> distances;
  std::optional<std::vector<double>> faceShear;
};

struct CommonSubdivisionDebugInfo {
  size_t vertexVertexCount = 0;
  size_t edgeTransverseCount = 0;
  size_t edgeParallelCount = 0;
  size_t faceVertexCount = 0;
  size_t edgeVertexCount = 0;
  size_t missingInputVertexCount = 0;
  size_t invalidPointsAlongAEndpointCount = 0;
  size_t invalidPointsAlongBEndpointCount = 0;
  size_t nonMonotonePointsAlongACount = 0;
  size_t nonMonotonePointsAlongBCount = 0;
  size_t emptyPointsAlongACount = 0;
  size_t emptyPointsAlongBCount = 0;
  size_t rawSubdivisionPointCount = 0;
  size_t expectedConstructedVertexCount = 0;
  size_t expectedConstructedEdgeCount = 0;
  size_t expectedConstructedFaceCount = 0;
  bool attemptedMeshConstruction = false;
  bool meshConstructed = false;
  size_t constructedVertexCount = 0;
  size_t constructedFaceCount = 0;
  std::optional<std::string> constructionError;
};

class GeoDrapeSolver {
public:
  explicit GeoDrapeSolver(SurfaceMeshData meshData,
                          const SignedHeatSolveOptions& heatOptions = {});
  GeoDrapeSolver(SurfaceMeshData meshData,
                 const SignedHeatSolveOptions& heatOptions,
                 const IntrinsicConstructionOptions& intrinsicOptions,
                 const RefinementOptions& refinementOptions = {});

  DrapeResult solve(const Vec2& seedXY,
                    double fabricAngle,
                    const DrapeSolveOptions& solveOptions = {});
  DrapeResult solve(const Vec2& seedXY,
                    double fabricAngle,
                    double fiberAngle,
                    const DrapeSolveOptions& solveOptions = {});
  DrapeResult solveFromIntrinsic(const SurfaceReference& seed,
                                 const TangentVectorRef& fabricDirection,
                                 double fiberAngle = 90.0,
                                 const DrapeSolveOptions& solveOptions = {});
  DrapeResult retrieve(ResultDomain retrieval = ResultDomain::Extrinsic,
                       bool sampleVertexShear = true);
  CommonSubdivisionDebugInfo debugCommonSubdivision(bool attemptMeshConstruction = false);

private:
  IntrinsicSolveInput adaptExtrinsicInput(const Vec2& seedXY,
                                          double fabricAngle,
                                          double fiberAngle,
                                          DrapeSolveMode mode,
                                          const TraceSettings& trace);
  CoreIntrinsicResult solveCore(const IntrinsicSolveInput& input);
  DrapeResult retrieveFromCore(const CoreIntrinsicResult& core,
                               ResultDomain retrieval,
                               bool sampleVertexShear);

  ReferenceGeometry reference_;
  ActiveIntrinsicDomain activeDomain_;
  TraceSettings traceDefaults_;
  SignedHeatSolveOptions heatOptions_;
  std::unique_ptr<CustomSignedHeatSolver> customHeatSolver_;
  std::optional<CoreIntrinsicResult> lastIntrinsicResult_;
  bool inputConnectivityPreserved_ = true;
};

DrapeResult solveDrape(const SurfaceMeshData& mesh,
                       const Vec2& seedXY,
                       double fabricAngle,
                       const SignedHeatSolveOptions& heatOptions = {},
                       const DrapeSolveOptions& solveOptions = {});

} // namespace geodesic_draping

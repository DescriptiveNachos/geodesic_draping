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
  bool sampleSecondaryShear = false;
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

  std::array<FaceHeatDirectionField, 2> faceDirections;
  std::optional<std::array<std::vector<double>, 2>> distances;
  std::optional<std::vector<double>> faceShearAnglesDegrees;
  std::optional<std::vector<double>> vertexShearAnglesDegrees;
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
  SeedProjection seed;
  std::array<Vec3, 4> directions;
  std::array<GeneratorTrace, 4> generators;
  SourceCurves sourceCurves;
  std::array<CustomSignedHeatResult, 2> customHeatSolves;
  std::array<FaceHeatDirectionField, 2> faceDirections;
  std::optional<std::array<std::vector<double>, 2>> distances;
  std::optional<std::array<std::vector<Vec3>, 2>> gradients;
  std::optional<std::vector<double>> faceShearAnglesDegrees;
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
                    double angleDegrees,
                    const DrapeSolveOptions& solveOptions = {});
  DrapeResult solve(const Vec2& seedXY,
                    double fabricAngleDegrees,
                    double fiberAngleDegrees,
                    const DrapeSolveOptions& solveOptions = {});
  DrapeResult retrieve(ResultDomain retrieval = ResultDomain::Extrinsic,
                       bool sampleVertexShear = true);

private:
  IntrinsicSolveInput adaptExtrinsicInput(const Vec2& seedXY,
                                          double fabricAngleDegrees,
                                          double fiberAngleDegrees,
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
  bool preservesInputConnectivity_ = true;
};

DrapeResult solveDrape(const SurfaceMeshData& mesh,
                       const Vec2& seedXY,
                       double angleDegrees,
                       const SignedHeatSolveOptions& heatOptions = {},
                       const DrapeSolveOptions& solveOptions = {});

} // namespace geodesic_draping

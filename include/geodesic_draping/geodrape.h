#pragma once

#include "geodesic_draping/custom_signed_heat.h"
#include "geodesic_draping/seed_projection.h"
#include "geodesic_draping/surface_construction.h"

#include "geometrycentral/surface/intrinsic_triangulation.h"
#include "geometrycentral/surface/barycentric_vector.h"
#include "geometrycentral/surface/vertex_position_geometry.h"
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

enum class RefinementMode {
  None,
  DelaunayFlip,
  DelaunayRefine,
};

enum class IntrinsicTriangulationBackend {
  Signpost,
  IntegerCoordinates,
};

struct IntrinsicConstructionOptions {
  IntrinsicTriangulationBackend backend = IntrinsicTriangulationBackend::Signpost;
};

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

struct TraceSettings {
  double traceLength = 10000.0;
  size_t maxIterations = geometrycentral::INVALID_IND;
};

struct DrapeSolveOptions {
  DrapeSolveMode mode = DrapeSolveMode::Complete;
  double fiberAngle = 90.0;
  AdvancedSolveOptions advanced;
};

struct IntrinsicGeneratorTrace {
  std::vector<geometrycentral::surface::SurfacePoint> points;
  bool hitBoundary = false;
  double length = 0.0;
};

struct IntrinsicSolveInput {
  geometrycentral::surface::SurfacePoint seed;
  std::array<geometrycentral::surface::BarycentricVector, 4> directions;
  DrapeSolveMode mode = DrapeSolveMode::Complete;
  TraceSettings trace;
};

struct CoreIntrinsicResult {
  DrapeSolveMode mode = DrapeSolveMode::Complete;
  geometrycentral::surface::SurfacePoint intrinsicSeed;
  std::array<geometrycentral::surface::BarycentricVector, 4> intrinsicDirections;
  std::array<IntrinsicGeneratorTrace, 4> generators;
  std::array<FaceHeatDirectionField, 2> directions;
  std::optional<std::array<std::vector<double>, 2>> distances;
  std::optional<std::vector<double>> faceShear;
};

struct ExtrinsicGeneratorTrace {
  std::vector<geometrycentral::Vector3> points;
  bool hitBoundary = false;
  double length = 0.0;
};

struct IntrinsicDrapeResult {
  const geometrycentral::surface::SurfaceMesh* mesh = nullptr;
  const geometrycentral::surface::IntrinsicGeometryInterface* geometry = nullptr;
  DrapeSolveMode mode = DrapeSolveMode::Complete;
  geometrycentral::surface::SurfacePoint seed;
  std::array<geometrycentral::surface::BarycentricVector, 4> directions;
  std::array<IntrinsicGeneratorTrace, 4> generators;
  std::array<geometrycentral::surface::FaceData<geometrycentral::Vector3>, 2> directionFields;
  std::optional<std::array<geometrycentral::surface::VertexData<double>, 2>> distances;
  std::optional<geometrycentral::surface::FaceData<double>> faceShear;
  std::optional<geometrycentral::surface::VertexData<double>> vertexShear;
};

struct ExtrinsicDrapeResult {
  const geometrycentral::surface::SurfaceMesh* mesh = nullptr;
  const geometrycentral::surface::VertexPositionGeometry* geometry = nullptr;
  DrapeSolveMode mode = DrapeSolveMode::Complete;
  geometrycentral::Vector3 seed = geometrycentral::Vector3::zero();
  std::array<geometrycentral::Vector3, 4> directions;
  std::array<ExtrinsicGeneratorTrace, 4> generators;
  std::optional<std::array<geometrycentral::surface::VertexData<double>, 2>> distances;
  std::optional<geometrycentral::surface::VertexData<double>> vertexShear;
};

struct SubdivisionDrapeResult {
  const geometrycentral::surface::SurfaceMesh* mesh = nullptr;
  DrapeSolveMode mode = DrapeSolveMode::Complete;
  geometrycentral::surface::VertexData<geometrycentral::Vector3> vertexPositions;
  geometrycentral::Vector3 seed = geometrycentral::Vector3::zero();
  std::array<geometrycentral::Vector3, 4> directions;
  std::array<ExtrinsicGeneratorTrace, 4> generators;
  std::array<geometrycentral::surface::FaceData<geometrycentral::Vector3>, 2> directionFields;
  std::optional<std::array<geometrycentral::surface::VertexData<double>, 2>> distances;
  std::optional<geometrycentral::surface::FaceData<double>> faceShear;
  std::optional<geometrycentral::surface::VertexData<double>> vertexShear;
};

class GeoDrapeSolver {
public:
  explicit GeoDrapeSolver(SurfaceMeshData meshData,
                          const SignedHeatSolveOptions& heatOptions = {});
  GeoDrapeSolver(SurfaceMeshData meshData,
                 const SignedHeatSolveOptions& heatOptions,
                 const IntrinsicConstructionOptions& intrinsicOptions,
                 const RefinementOptions& refinementOptions = {});

  const CoreIntrinsicResult& solve(const Vec2& seedXY,
                                   double fabricAngle,
                                   const DrapeSolveOptions& solveOptions = {});
  const CoreIntrinsicResult& solveFromIntrinsic(
      const geometrycentral::surface::SurfacePoint& seed,
      const geometrycentral::surface::BarycentricVector& fabricDirection,
      double fiberAngle = 90.0,
      const DrapeSolveOptions& solveOptions = {});
  const CoreIntrinsicResult& lastResult() const;
  IntrinsicDrapeResult retrieveIntrinsic(bool sampleVertexShear = false) const;
  ExtrinsicDrapeResult retrieveExtrinsic(bool sampleVertexShear = false) const;
  SubdivisionDrapeResult retrieveSubdivision(bool sampleVertexShear = false) const;

private:
  IntrinsicSolveInput adaptExtrinsicInput(const Vec2& seedXY,
                                          double fabricAngle,
                                          double fiberAngle,
                                          DrapeSolveMode mode,
                                          const TraceSettings& trace);
  CoreIntrinsicResult solveCore(const IntrinsicSolveInput& input);

  SurfaceMeshData meshData_;
  GeometryCentralSurface inputSurface_;
  std::unique_ptr<geometrycentral::surface::IntrinsicTriangulation> intrinsicTriangulation_;
  bool useCommonSubdivisionInputAdapter_ = false;
  TraceSettings traceDefaults_;
  SignedHeatSolveOptions heatOptions_;
  std::unique_ptr<CustomSignedHeatSolver> customHeatSolver_;
  std::optional<CoreIntrinsicResult> lastIntrinsicResult_;
  bool inputConnectivityPreserved_ = true;
};

} // namespace geodesic_draping

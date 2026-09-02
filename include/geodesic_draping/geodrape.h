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

enum class RetrievalDomain {
  Intrinsic,
  Extrinsic,
  Subdivision,
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

struct TraceSettings {
  double traceLength = 10000.0;
  size_t maxIterations = geometrycentral::INVALID_IND;
};

struct DrapeSolveOptions {
  DrapeSolveMode mode = DrapeSolveMode::Complete;
  double fiberAngle = 90.0;
  AdvancedTraceOptions trace;
};

struct RetrievalOptions {
  RetrievalDomain domain = RetrievalDomain::Extrinsic;
  bool sampleVertexShear = false;
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

struct DrapeResult {
  RetrievalDomain domain = RetrievalDomain::Extrinsic;
  DrapeSolveMode mode = DrapeSolveMode::Complete;

  // Keeps solver-owned geometry alive for one-shot solveDrape() results.
  std::shared_ptr<const void> storageOwner;

  const geometrycentral::surface::SurfaceMesh* mesh = nullptr;
  const geometrycentral::surface::IntrinsicGeometryInterface* intrinsicGeometry = nullptr;
  const geometrycentral::surface::VertexPositionGeometry* extrinsicGeometry = nullptr;

  std::optional<geometrycentral::surface::VertexData<geometrycentral::Vector3>> vertexPositions;

  std::optional<geometrycentral::surface::SurfacePoint> intrinsicSeed;
  std::optional<std::array<geometrycentral::surface::BarycentricVector, 4>> intrinsicDirections;
  std::optional<std::array<IntrinsicGeneratorTrace, 4>> intrinsicGenerators;

  std::optional<geometrycentral::Vector3> extrinsicSeed;
  std::optional<std::array<geometrycentral::Vector3, 4>> extrinsicDirections;
  std::optional<std::array<ExtrinsicGeneratorTrace, 4>> extrinsicGenerators;

  std::optional<std::array<geometrycentral::surface::FaceData<geometrycentral::Vector3>, 2>> directionFields;
  std::optional<std::array<geometrycentral::surface::VertexData<double>, 2>> distances;
  std::optional<geometrycentral::surface::FaceData<double>> faceShear;
  std::optional<geometrycentral::surface::VertexData<double>> vertexShear;
};

class GeoDrapeSolver {
public:
  explicit GeoDrapeSolver(
      SurfaceMeshData meshData,
      const SignedHeatSolveOptions& heatOptions = {},
      IntrinsicTriangulationBackend intrinsicBackend = IntrinsicTriangulationBackend::Signpost,
      const RefinementOptions& refinementOptions = {});

  DrapeResult solve(const Vec2& seedXY,
                    double fabricAngle,
                    const DrapeSolveOptions& solveOptions = {},
                    const RetrievalOptions& retrievalOptions = {});
  DrapeResult solveFromIntrinsic(
      const geometrycentral::surface::SurfacePoint& seed,
      const geometrycentral::surface::BarycentricVector& fabricDirection,
      const DrapeSolveOptions& solveOptions = {},
      const RetrievalOptions& retrievalOptions = {});
  DrapeResult retrieve(const RetrievalOptions& retrievalOptions = {}) const;
  const CoreIntrinsicResult& lastCoreResult() const;

private:
  IntrinsicSolveInput adaptExtrinsicInput(const Vec2& seedXY,
                                          double fabricAngle,
                                          double fiberAngle,
                                          DrapeSolveMode mode,
                                          const TraceSettings& trace);
  IntrinsicSolveInput adaptIntrinsicInput(
      const geometrycentral::surface::SurfacePoint& seed,
      const geometrycentral::surface::BarycentricVector& fabricDirection,
      const DrapeSolveOptions& solveOptions,
      const TraceSettings& trace);
  CoreIntrinsicResult solveCore(const IntrinsicSolveInput& input);
  DrapeResult retrieveIntrinsic(bool sampleVertexShear = false) const;
  DrapeResult retrieveExtrinsic(bool sampleVertexShear = false) const;
  DrapeResult retrieveSubdivision(bool sampleVertexShear = false) const;

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

DrapeResult solveDrape(SurfaceMeshData meshData,
                       const Vec2& seedXY,
                       double fabricAngle,
                       const SignedHeatSolveOptions& heatOptions = {},
                       const DrapeSolveOptions& solveOptions = {},
                       const RetrievalOptions& retrievalOptions = {},
                       IntrinsicTriangulationBackend intrinsicBackend = IntrinsicTriangulationBackend::Signpost,
                       const RefinementOptions& refinementOptions = {});

} // namespace geodesic_draping

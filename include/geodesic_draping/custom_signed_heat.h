#pragma once

#include "geodesic_draping/signed_heat.h"

#include <array>
#include <complex>
#include <memory>
#include <vector>

namespace geodesic_draping {

using EdgeHeatField = std::vector<std::complex<double>>;
using FaceHeatDirectionField = std::vector<geometrycentral::Vector3>;

enum class VertexDirectionAveraging {
  FaceArea,
  CornerAngle,
};

struct DiffusedHeatFieldResult {
  std::vector<geometrycentral::surface::Curve> preprocessedSourceCurves;
  EdgeHeatField sourceEdgeHeatField;
  EdgeHeatField diffusedEdgeHeatField;
};

struct CustomSignedHeatResult {
  DiffusedHeatFieldResult diffusion;
  FaceHeatDirectionField normalizedFaceDirections;
  std::vector<double> distance;
};

struct CustomSignedHeatStageTimings {
  double curveConversionSeconds = 0.0;
  double preprocessSeconds = 0.0;
  double sourceSeconds = 0.0;
  double diffuseSeconds = 0.0;
  double normalizeSeconds = 0.0;
  double distanceSeconds = 0.0;
  double totalSeconds = 0.0;
  bool heatSolverAlreadyInitialized = false;
  bool poissonSolverAlreadyInitialized = false;
};

struct TimedCustomSignedHeatResult {
  CustomSignedHeatResult result;
  CustomSignedHeatStageTimings timings;
};

class CustomSignedHeatSolver {
public:
  CustomSignedHeatSolver(GeometryCentralSurface& surface, double diffusionTimeCoefficient = 1.0);
  CustomSignedHeatSolver(geometrycentral::surface::SurfaceMesh& mesh,
                         geometrycentral::surface::IntrinsicGeometryInterface& geometry,
                         double diffusionTimeCoefficient = 1.0);

  DiffusedHeatFieldResult solveDiffusedEdgeHeatField(
      const geometrycentral::surface::Curve& sourceCurve,
      const SignedHeatSolveOptions& options = {});
  DiffusedHeatFieldResult solveDiffusedEdgeHeatField(
      const std::vector<SurfaceReference>& sourceCurve,
      const SignedHeatSolveOptions& options = {});

  CustomSignedHeatResult solve(const geometrycentral::surface::Curve& sourceCurve,
                               const SignedHeatSolveOptions& options = {},
                               bool computeDistance = false);
  CustomSignedHeatResult solve(const std::vector<SurfaceReference>& sourceCurve,
                               const SignedHeatSolveOptions& options = {},
                               bool computeDistance = false);
  TimedCustomSignedHeatResult solveTimed(const geometrycentral::surface::Curve& sourceCurve,
                                         const SignedHeatSolveOptions& options = {},
                                         bool computeDistance = false);
  TimedCustomSignedHeatResult solveTimed(const std::vector<SurfaceReference>& sourceCurve,
                                         const SignedHeatSolveOptions& options = {},
                                         bool computeDistance = false);

  std::array<CustomSignedHeatResult, 2> solve(
      const std::array<geometrycentral::surface::Curve, 2>& sourceCurves,
      const SignedHeatSolveOptions& options = {},
      bool computeDistance = false);
  std::array<CustomSignedHeatResult, 2> solve(const SourceCurves& sourceCurves,
                                              const SignedHeatSolveOptions& options = {},
                                              bool computeDistance = false);
  std::array<TimedCustomSignedHeatResult, 2> solveTimed(
      const std::array<geometrycentral::surface::Curve, 2>& sourceCurves,
      const SignedHeatSolveOptions& options = {},
      bool computeDistance = false);
  std::array<TimedCustomSignedHeatResult, 2> solveTimed(const SourceCurves& sourceCurves,
                                                        const SignedHeatSolveOptions& options = {},
                                                        bool computeDistance = false);

private:
  geometrycentral::surface::SurfaceMesh& mesh_;
  geometrycentral::surface::IntrinsicGeometryInterface& geometry_;

  double shortTime_ = 0.0;
  double meanNodeDistance_ = 0.0;

  geometrycentral::SparseMatrix<double> massMatrix_;
  geometrycentral::SparseMatrix<double> doubleMassMatrix_;
  geometrycentral::SparseMatrix<double> doubleConnectionLaplacian_;
  geometrycentral::SparseMatrix<double> doubleVectorOperator_;
  std::unique_ptr<geometrycentral::LinearSolver<std::complex<double>>> heatFieldSolver_;
  std::unique_ptr<geometrycentral::PositiveDefiniteSolver<double>> poissonSolver_;

  std::vector<geometrycentral::surface::Curve> preprocessCurves(
      const std::vector<geometrycentral::surface::Curve>& curves) const;
  EdgeHeatField buildSourceEdgeHeatField(
      const std::vector<geometrycentral::surface::Curve>& curves) const;
  EdgeHeatField diffuseEdgeHeatField(
      const EdgeHeatField& sourceEdgeHeatField,
      const std::vector<geometrycentral::surface::Curve>& curves,
      const SignedHeatSolveOptions& options);
  FaceHeatDirectionField sampleAndNormalizeFaceDirections(const EdgeHeatField& diffusedEdgeHeatField);
  std::vector<double> integrateVectorFieldToDistance(
      const FaceHeatDirectionField& normalizedFaceDirections,
      const std::vector<geometrycentral::surface::Curve>& curves,
      const SignedHeatSolveOptions& options);

  void ensureHaveHeatFieldSolver();
  void ensureHavePoissonSolver();
  geometrycentral::SparseMatrix<double> buildCrouzeixRaviartDoubleConnectionLaplacian() const;
  geometrycentral::SparseMatrix<double> buildCrouzeixRaviartDoubleMassMatrix() const;

  void buildSignedCurveSource(const geometrycentral::surface::Curve& curve,
                              EdgeHeatField& sourceEdgeHeatField) const;
  double lengthOfSegment(const geometrycentral::surface::SurfacePoint& pA,
                         const geometrycentral::surface::SurfacePoint& pB) const;
  geometrycentral::surface::SurfacePoint midSegmentSurfacePoint(
      const geometrycentral::surface::SurfacePoint& pA,
      const geometrycentral::surface::SurfacePoint& pB) const;
  std::complex<double> projectedNormal(const geometrycentral::surface::SurfacePoint& pA,
                                       const geometrycentral::surface::SurfacePoint& pB,
                                       const geometrycentral::surface::Edge& edge) const;
  double scalarCrouzeixRaviart(const geometrycentral::surface::SurfacePoint& point,
                               const geometrycentral::surface::Edge& edge) const;
  double averageValueOnSource(const geometrycentral::Vector<double>& phi,
                              const std::vector<geometrycentral::surface::Curve>& curves) const;
};

FaceHeatDirectionField sampleAndNormalizeFaceDirections(GeometryCentralSurface& surface,
                                                        const EdgeHeatField& diffusedEdgeHeatField);

std::vector<Vec3> directionsToExtrinsicVectors(
    GeometryCentralSurface& surface,
    const FaceHeatDirectionField& normalizedFaceDirections);

std::vector<double> computeFaceShearAnglesDegrees(
    GeometryCentralSurface& surface,
    const FaceHeatDirectionField& normalizedFaceDirections0,
    const FaceHeatDirectionField& normalizedFaceDirections1);
std::vector<double> computeFaceShearAnglesDegrees(
    geometrycentral::surface::SurfaceMesh& mesh,
    geometrycentral::surface::IntrinsicGeometryInterface& geometry,
    const FaceHeatDirectionField& normalizedFaceDirections0,
    const FaceHeatDirectionField& normalizedFaceDirections1);

CustomSignedHeatResult computeCustomSignedHeat(GeometryCentralSurface& surface,
                                                         const std::vector<SurfaceReference>& sourceCurve,
                                                         const SignedHeatSolveOptions& options = {});

std::array<CustomSignedHeatResult, 2> computeCustomSignedHeat(GeometryCentralSurface& surface,
                                                                        const SourceCurves& sourceCurves,
                                                                        const SignedHeatSolveOptions& options = {});

} // namespace geodesic_draping

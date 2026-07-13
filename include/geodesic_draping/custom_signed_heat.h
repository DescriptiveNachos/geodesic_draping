#pragma once

#include "geometrycentral/surface/signed_heat_method.h"

#include <array>
#include <complex>
#include <memory>
#include <vector>

namespace geodesic_draping {

using EdgeHeatField = std::vector<std::complex<double>>;
using FaceHeatDirectionField = std::vector<geometrycentral::Vector3>;

struct SignedHeatSolveOptions {
  bool preserveSourceNormals = false;
  geometrycentral::LevelSetConstraint levelSetConstraint = geometrycentral::LevelSetConstraint::None;
  double softLevelSetWeight = -1.0;
  double diffusionTimeCoefficient = 1.0;
};

struct DiffusedHeatFieldResult {
  std::vector<geometrycentral::surface::Curve> preprocessedCurves;
  EdgeHeatField sourceEdgeHeatField;
  EdgeHeatField diffusedEdgeHeatField;
};

struct CustomSignedHeatResult {
  DiffusedHeatFieldResult diffusion;
  FaceHeatDirectionField normalizedFaceDirections;
  std::vector<double> distance;
};

class CustomSignedHeatSolver {
public:
  CustomSignedHeatSolver(geometrycentral::surface::SurfaceMesh& mesh,
                         geometrycentral::surface::IntrinsicGeometryInterface& geometry,
                         double diffusionTimeCoefficient = 1.0);

  DiffusedHeatFieldResult solveDiffusedEdgeHeatField(
      const geometrycentral::surface::Curve& sourceCurve,
      const SignedHeatSolveOptions& options = {});

  CustomSignedHeatResult solve(const geometrycentral::surface::Curve& sourceCurve,
                               const SignedHeatSolveOptions& options = {},
                               bool computeDistance = false);

  std::array<CustomSignedHeatResult, 2> solve(
      const std::array<geometrycentral::surface::Curve, 2>& sourceCurves,
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

std::vector<double> computeFaceShearAnglesDegrees(
    geometrycentral::surface::SurfaceMesh& mesh,
    geometrycentral::surface::IntrinsicGeometryInterface& geometry,
    const FaceHeatDirectionField& normalizedFaceDirections0,
    const FaceHeatDirectionField& normalizedFaceDirections1);

} // namespace geodesic_draping

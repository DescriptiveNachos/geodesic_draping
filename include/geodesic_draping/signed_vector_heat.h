#pragma once

#include "geodesic_draping/signed_heat.h"

#include <array>
#include <complex>
#include <memory>
#include <vector>

namespace geodesic_draping {

using EdgeVectorHeat = std::vector<std::complex<double>>;
using FaceVectorHeat = std::vector<Vec3>;

struct DiffusedVectorHeatResult {
  std::vector<std::vector<SurfaceReference>> preprocessedSourceCurves;
  EdgeVectorHeat sourceEdgeVectorHeat;
  EdgeVectorHeat diffusedEdgeVectorHeat;
};

struct SignedVectorHeatResult {
  DiffusedVectorHeatResult diffusion;
  FaceVectorHeat normalizedFaceVectorHeat;
  std::vector<Vec3> vertexVectorHeat;
};

class SignedVectorHeatSolver {
public:
  SignedVectorHeatSolver(GeometryCentralSurface& surface, double diffusionTimeCoefficient = 1.0);

  DiffusedVectorHeatResult solveDiffusedEdgeVectorHeat(
      const std::vector<SurfaceReference>& sourceCurve,
      const SignedHeatSolveOptions& options = {});

private:
  geometrycentral::surface::SurfaceMesh& mesh_;
  geometrycentral::surface::IntrinsicGeometryInterface& geometry_;

  double shortTime_ = 0.0;
  double meanNodeDistance_ = 0.0;

  geometrycentral::SparseMatrix<double> massMatrix_;
  geometrycentral::SparseMatrix<double> doubleMassMatrix_;
  geometrycentral::SparseMatrix<double> doubleConnectionLaplacian_;
  geometrycentral::SparseMatrix<double> doubleVectorOperator_;
  std::unique_ptr<geometrycentral::LinearSolver<std::complex<double>>> vectorHeatSolver_;

  std::vector<geometrycentral::surface::Curve> preprocessCurves(
      const std::vector<geometrycentral::surface::Curve>& curves) const;
  EdgeVectorHeat buildSourceEdgeVectorHeat(
      const std::vector<geometrycentral::surface::Curve>& curves) const;
  EdgeVectorHeat diffuseEdgeVectorHeat(
      const EdgeVectorHeat& sourceEdgeVectorHeat,
      const std::vector<geometrycentral::surface::Curve>& curves,
      const SignedHeatSolveOptions& options);

  void ensureHaveVectorHeatSolver();
  geometrycentral::SparseMatrix<double> buildCrouzeixRaviartDoubleConnectionLaplacian() const;
  geometrycentral::SparseMatrix<double> buildCrouzeixRaviartDoubleMassMatrix() const;

  void buildSignedCurveSource(const geometrycentral::surface::Curve& curve,
                              EdgeVectorHeat& sourceEdgeVectorHeat) const;
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
};

FaceVectorHeat sampleAndNormalizeFaceVectorHeat(GeometryCentralSurface& surface,
                                                const EdgeVectorHeat& diffusedEdgeVectorHeat);

std::vector<Vec3> averageNormalizedFaceVectorHeatToVerticesReference(
    GeometryCentralSurface& surface,
    const FaceVectorHeat& normalizedFaceVectorHeat);

SignedVectorHeatResult computeSignedVectorHeat(GeometryCentralSurface& surface,
                                               const std::vector<SurfaceReference>& sourceCurve,
                                               const SignedHeatSolveOptions& options = {});

std::array<SignedVectorHeatResult, 2> computeSignedVectorHeats(GeometryCentralSurface& surface,
                                                               const SourceCurves& sourceCurves,
                                                               const SignedHeatSolveOptions& options = {});

} // namespace geodesic_draping

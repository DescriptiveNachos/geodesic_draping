#include "geodesic_draping/signed_heat.h"

#include <stdexcept>

namespace geodesic_draping {
namespace {

geometrycentral::SignedHeatOptions toGeometryCentralOptions(const SignedHeatSolveOptions& options) {
  geometrycentral::SignedHeatOptions gcOptions;
  gcOptions.preserveSourceNormals = options.preserveSourceNormals;
  gcOptions.levelSetConstraint = options.levelSetConstraint;
  gcOptions.softLevelSetWeight = options.softLevelSetWeight;
  return gcOptions;
}

std::vector<double> toOrderedVector(const geometrycentral::surface::VertexData<double>& data,
                                    geometrycentral::surface::SurfaceMesh& mesh) {
  std::vector<double> values(mesh.nVertices(), 0.0);
  for (geometrycentral::surface::Vertex v : mesh.vertices()) {
    values[v.getIndex()] = data[v];
  }
  return values;
}

geometrycentral::surface::SurfaceMesh& requireMesh(GeometryCentralSurface& surface) {
  if (!surface.mesh) {
    throw std::runtime_error("SignedHeatDistanceSolver requires a valid geometry-central mesh");
  }
  return *surface.mesh;
}

geometrycentral::surface::IntrinsicGeometryInterface& requireGeometry(GeometryCentralSurface& surface) {
  if (!surface.geometry) {
    throw std::runtime_error("SignedHeatDistanceSolver requires a valid geometry-central geometry");
  }
  return *surface.geometry;
}

} // namespace

geometrycentral::surface::SurfacePoint toGeometryCentralSurfacePoint(
    geometrycentral::surface::SurfaceMesh& mesh,
    const SurfaceReference& ref) {
  switch (ref.type) {
  case SurfaceReferenceType::Vertex:
    return geometrycentral::surface::SurfacePoint(mesh.vertex(ref.elementIndex));
  case SurfaceReferenceType::Edge:
    if (ref.params.size() < 1) {
      throw std::runtime_error("edge surface reference requires one parameter");
    }
    return geometrycentral::surface::SurfacePoint(mesh.edge(ref.elementIndex), ref.params[0]);
  case SurfaceReferenceType::Face: {
    if (ref.params.size() < 2) {
      throw std::runtime_error("face surface reference requires at least two barycentric parameters");
    }
    const double z = ref.params.size() >= 3 ? ref.params[2] : 1.0 - ref.params[0] - ref.params[1];
    return geometrycentral::surface::SurfacePoint(
        mesh.face(ref.elementIndex),
        geometrycentral::Vector3{ref.params[0], ref.params[1], z});
  }
  }
  throw std::runtime_error("unknown surface reference type");
}

geometrycentral::surface::Curve toGeometryCentralCurve(
    geometrycentral::surface::SurfaceMesh& mesh,
    const std::vector<SurfaceReference>& refs,
    bool isSigned) {
  geometrycentral::surface::Curve curve;
  curve.isSigned = isSigned;
  curve.nodes.reserve(refs.size());
  for (const SurfaceReference& ref : refs) {
    curve.nodes.push_back(toGeometryCentralSurfacePoint(mesh, ref));
  }
  return curve;
}

std::vector<double> computeSignedHeatDistance(GeometryCentralSurface& surface,
                                              const std::vector<SurfaceReference>& sourceCurve,
                                              const SignedHeatSolveOptions& options) {
  if (!surface.mesh || !surface.geometry) {
    throw std::runtime_error("computeSignedHeatDistance requires a valid geometry-central surface");
  }

  SignedHeatDistanceSolver solver(surface, options.diffusionTimeCoefficient);
  return solver.computeDistance(sourceCurve, options);
}

std::array<std::vector<double>, 2> computeSignedHeatDistances(GeometryCentralSurface& surface,
                                                              const SourceCurves& sourceCurves,
                                                              const SignedHeatSolveOptions& options) {
  SignedHeatDistanceSolver solver(surface, options.diffusionTimeCoefficient);
  return solver.computeDistances(sourceCurves, options);
}

SignedHeatDistanceSolver::SignedHeatDistanceSolver(GeometryCentralSurface& surface,
                                                   double diffusionTimeCoefficient)
    : mesh_(requireMesh(surface)),
      solver_(requireGeometry(surface), diffusionTimeCoefficient) {}

std::vector<double> SignedHeatDistanceSolver::computeDistance(
    const std::vector<SurfaceReference>& sourceCurve,
    const SignedHeatSolveOptions& options) {
  std::vector<geometrycentral::surface::Curve> curves{
      toGeometryCentralCurve(mesh_, sourceCurve, true),
  };
  const geometrycentral::surface::VertexData<double> distance =
      solver_.computeDistance(curves, toGeometryCentralOptions(options));
  return toOrderedVector(distance, mesh_);
}

std::array<std::vector<double>, 2> SignedHeatDistanceSolver::computeDistances(
    const SourceCurves& sourceCurves,
    const SignedHeatSolveOptions& options) {
  return {
      computeDistance(sourceCurves.curves[0], options),
      computeDistance(sourceCurves.curves[1], options),
  };
}

} // namespace geodesic_draping

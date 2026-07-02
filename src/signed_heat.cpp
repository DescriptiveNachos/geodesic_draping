#include "geodesic_draping/signed_heat.h"

#include <stdexcept>

namespace geodesic_draping {

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

} // namespace geodesic_draping

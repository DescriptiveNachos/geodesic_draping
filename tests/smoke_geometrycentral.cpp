#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/trace_geodesic.h"
#include "geometrycentral/surface/vertex_position_geometry.h"
#include "geometrycentral/utilities/vector3.h"

#include <Eigen/Dense>

#include <cassert>
#include <cmath>

namespace {

bool near(double a, double b, double tolerance = 1e-10) {
  return std::abs(a - b) <= tolerance;
}

} // namespace

int main() {
  using geometrycentral::Vector2;
  using geometrycentral::Vector3;
  using geometrycentral::surface::Face;
  using geometrycentral::surface::ManifoldSurfaceMesh;
  using geometrycentral::surface::SurfacePoint;
  using geometrycentral::surface::TraceOptions;
  using geometrycentral::surface::VertexPositionGeometry;
  using geometrycentral::surface::traceGeodesic;

  Eigen::Matrix<int, 2, 3> faces;
  faces << 0, 1, 2,
           0, 2, 3;

  Eigen::Matrix<double, 4, 3> vertices;
  vertices << 0.0, 0.0, 0.0,
              1.0, 0.0, 0.0,
              1.0, 1.0, 0.0,
              0.0, 1.0, 0.0;

  ManifoldSurfaceMesh mesh(faces);
  VertexPositionGeometry geometry(mesh, vertices);

  geometry.requireFaceAreas();
  double totalArea = 0.0;
  for (Face f : mesh.faces()) {
    totalArea += geometry.faceAreas[f];
  }
  assert(near(totalArea, 1.0));
  geometry.unrequireFaceAreas();

  geometry.requireEdgeLengths();
  assert(mesh.nEdges() == 5);
  for (auto e : mesh.edges()) {
    assert(geometry.edgeLengths[e] > 0.0);
  }
  geometry.unrequireEdgeLengths();

  TraceOptions options;
  options.includePath = true;
  options.errorOnProblem = true;
  SurfacePoint start(mesh.face(0), Vector3{0.25, 0.25, 0.5});
  auto result = traceGeodesic(geometry, start, Vector2{0.1, 0.0}, options);

  assert(result.hasPath);
  assert(!result.pathPoints.empty());
  assert(result.length > 0.0);

  return 0;
}

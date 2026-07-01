#include "fixture_io.h"
#include "geodesic_draping/geometrycentral_adapter.h"

#include <cassert>
#include <cmath>
#include <filesystem>
#include <iostream>
#include <string>

namespace {

bool near(double a, double b, double tolerance = 1e-9) {
  return std::abs(a - b) <= tolerance;
}

void testFixture(const std::filesystem::path& root, const std::string& name) {
  const std::filesystem::path fixtureDir = root / name;
  const geodesic_draping::SurfaceMeshData meshData = geodesic_draping::fixture_io::loadMesh(fixtureDir);
  auto surface = geodesic_draping::makeGeometryCentralSurface(meshData);

  assert(surface.mesh != nullptr);
  assert(surface.geometry != nullptr);
  assert(surface.mesh->nVertices() == meshData.vertices.size());
  assert(surface.mesh->nFaces() == meshData.faces.size());

  surface.geometry->requireVertexPositions();
  for (geometrycentral::surface::Vertex v : surface.mesh->vertices()) {
    const auto p = surface.geometry->vertexPositions[v];
    const auto expected = meshData.vertices[v.getIndex()];
    assert(near(p.x, expected.x()));
    assert(near(p.y, expected.y()));
    assert(near(p.z, expected.z()));
  }
  surface.geometry->unrequireVertexPositions();

  surface.geometry->requireFaceAreas();
  double totalArea = 0.0;
  for (geometrycentral::surface::Face f : surface.mesh->faces()) {
    assert(surface.geometry->faceAreas[f] > 0.0);
    totalArea += surface.geometry->faceAreas[f];
  }
  assert(totalArea > 0.0);
  surface.geometry->unrequireFaceAreas();

  surface.geometry->requireEdgeLengths();
  for (geometrycentral::surface::Edge e : surface.mesh->edges()) {
    assert(surface.geometry->edgeLengths[e] > 0.0);
  }
  surface.geometry->unrequireEdgeLengths();
}

} // namespace

int main() {
  const std::filesystem::path fixtureRoot = GEODESIC_DRAPING_TEST_DATA_DIR;
  testFixture(fixtureRoot, "tiny_planar");
  testFixture(fixtureRoot, "small_curved");
  testFixture(fixtureRoot, "demo_part");
  return 0;
}

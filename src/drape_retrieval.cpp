#include "geodrape_internal.h"

#include <stdexcept>

namespace geodesic_draping {
namespace {

std::vector<Face> meshFaces(gcs::SurfaceMesh& mesh) {
  std::vector<Face> faces;
  faces.reserve(mesh.nFaces());
  for (gcs::Face face : mesh.faces()) {
    Face out{};
    size_t i = 0;
    for (gcs::Vertex vertex : face.adjacentVertices()) {
      out[i] = vertex.getIndex();
      ++i;
    }
    faces.push_back(out);
  }
  return faces;
}

size_t localSideIndex(gcs::Halfedge target) {
  size_t side = 0;
  for (gcs::Halfedge halfedge : target.face().adjacentHalfedges()) {
    if (halfedge == target) {
      return side;
    }
    ++side;
  }
  throw std::runtime_error("halfedge is not incident on its face");
}

std::vector<Vec3> toVec3Vector(const gcs::VertexData<geometrycentral::Vector3>& values) {
  std::vector<Vec3> out(values.getMesh()->nVertices(), Vec3::Zero());
  for (gcs::Vertex vertex : values.getMesh()->vertices()) {
    const geometrycentral::Vector3& value = values[vertex];
    out[vertex.getIndex()] = Vec3(value.x, value.y, value.z);
  }
  return out;
}

DrapeTrace makeTrace(const GeneratorTrace& generator, ResultDomain domain) {
  DrapeTrace trace;
  trace.hitBoundary = generator.hitBoundary;
  trace.length = generator.length;
  if (domain == ResultDomain::Intrinsic) {
    trace.intrinsicPoints = generator.surfaceReferences;
  } else {
    trace.extrinsicPoints = generator.points;
  }
  return trace;
}

} // namespace

SeedProjection toExtrinsicSeed(ReferenceGeometry& reference,
                               ActiveIntrinsicDomain& activeDomain,
                               const SurfaceReference& intrinsicSeed) {
  const gcs::SurfacePoint intrinsicPoint =
      toGeometryCentralSurfacePoint(activeDomain.mesh(), intrinsicSeed);
  const gcs::SurfacePoint inputPoint =
      activeDomain.intrinsicToInput(intrinsicPoint).inSomeFace();

  SeedProjection seed;
  seed.surfacePoint.faceIndex = inputPoint.face.getIndex();
  seed.surfacePoint.barycentric =
      Vec3(inputPoint.faceCoords.x, inputPoint.faceCoords.y, inputPoint.faceCoords.z);
  seed.cartesian = interpolateSurfacePoint(inputPoint, *reference.surface().geometry);
  return seed;
}

ResultMesh makeExtrinsicResultMesh(const SurfaceMeshData& meshData) {
  ResultMesh mesh;
  mesh.domain = ResultDomain::Extrinsic;
  mesh.faces = meshData.faces;
  mesh.vertices3D = meshData.vertices;
  return mesh;
}

ResultMesh makeIntrinsicResultMesh(gcs::SurfaceMesh& mesh,
                                   gcs::IntrinsicGeometryInterface& geometry) {
  ResultMesh out;
  out.domain = ResultDomain::Intrinsic;
  out.faces = meshFaces(mesh);
  out.edgeLengths = std::vector<std::array<double, 3>>{};
  out.gluingMap = std::vector<FaceGluingMap>{};
  out.edgeLengths->reserve(mesh.nFaces());
  out.gluingMap->reserve(mesh.nFaces());

  geometry.requireEdgeLengths();
  for (gcs::Face face : mesh.faces()) {
    std::array<double, 3> lengths{};
    FaceGluingMap gluing{};
    size_t side = 0;
    for (gcs::Halfedge halfedge : face.adjacentHalfedges()) {
      lengths[side] = geometry.edgeLengths[halfedge.edge()];
      if (halfedge.twin().isInterior()) {
        const gcs::Halfedge twin = halfedge.twin();
        gluing[side] = {
            static_cast<int>(twin.face().getIndex()),
            static_cast<int>(localSideIndex(twin)),
        };
      } else {
        gluing[side] = {-1, -1};
      }
      ++side;
    }
    out.edgeLengths->push_back(lengths);
    out.gluingMap->push_back(gluing);
  }
  geometry.unrequireEdgeLengths();
  return out;
}

ResultMesh makeSubdivisionResultMesh(gcs::CommonSubdivision& subdivision,
                                     gcs::VertexPositionGeometry& inputGeometry) {
  subdivision.constructMesh();
  ResultMesh out;
  out.domain = ResultDomain::Subdivision;
  out.faces = meshFaces(*subdivision.mesh);
  inputGeometry.requireVertexPositions();
  out.vertices3D = toVec3Vector(subdivision.interpolateAcrossA(inputGeometry.vertexPositions));
  inputGeometry.unrequireVertexPositions();
  return out;
}

std::array<TraceFamily, 2> makeTraceFamilies(const std::array<GeneratorTrace, 4>& generators,
                                             ResultDomain domain) {
  return {
      TraceFamily{makeTrace(generators[0], domain), makeTrace(generators[1], domain)},
      TraceFamily{makeTrace(generators[2], domain), makeTrace(generators[3], domain)},
  };
}

DrapeResult GeoDrapeSolver::retrieveFromCore(const CoreIntrinsicResult& core,
                                             ResultDomain retrieval,
                                             bool sampleVertexShear) {
  DrapeResult result;
  result.domain = retrieval;
  result.mode = core.mode;
  const bool needsExtrinsicOrigin =
      retrieval == ResultDomain::Extrinsic || retrieval == ResultDomain::Subdivision;
  const std::optional<SeedProjection> extrinsicSeed =
      needsExtrinsicOrigin
          ? std::optional<SeedProjection>{toExtrinsicSeed(reference_, activeDomain_, core.intrinsicSeed)}
          : std::nullopt;

  if (retrieval == ResultDomain::Intrinsic) {
    result.mesh = makeIntrinsicResultMesh(activeDomain_.mesh(), activeDomain_.geometry());
    result.origin.intrinsicPoint = core.intrinsicSeed;
    result.origin.intrinsicFamilyDirections = {core.intrinsicDirections[0], core.intrinsicDirections[2]};
    result.traces = makeTraceFamilies(core.generators, ResultDomain::Intrinsic);
    result.directions = core.directions;
    result.faceShear = core.faceShear;
    result.distances = core.distances;
    if (sampleVertexShear && core.faceShear) {
      result.vertexShear =
          averageIntrinsicFaceScalarsToVertices(activeDomain_.mesh(), *core.faceShear);
    }
  } else if (retrieval == ResultDomain::Subdivision) {
    gcs::CommonSubdivision& subdivision = activeDomain_.triangulation().getCommonSubdivision();
    subdivision.constructMesh();
    result.mesh = makeSubdivisionResultMesh(subdivision, *reference_.surface().geometry);
    result.mesh.intrinsicSourceFaceColor =
        faceDataToVector(subdivision.copyFromB(gcs::niceColors(activeDomain_.mesh())));
    result.origin.intrinsicPoint = core.intrinsicSeed;
    result.origin.intrinsicFamilyDirections = {core.intrinsicDirections[0], core.intrinsicDirections[2]};
    result.origin.extrinsicPoint = extrinsicSeed->cartesian;
    result.origin.extrinsicFamilyDirections = {core.cartesianDirections[0], core.cartesianDirections[2]};
    result.traces = makeTraceFamilies(core.generators, ResultDomain::Subdivision);

    result.directions = {
        faceVectorDataToVector(subdivision.copyFromB(
            activeFaceVectorData(activeDomain_.mesh(), core.directions[0]))),
        faceVectorDataToVector(subdivision.copyFromB(
            activeFaceVectorData(activeDomain_.mesh(), core.directions[1]))),
    };
    if (core.faceShear) {
      result.faceShear = faceDataToVector(subdivision.copyFromB(
          activeFaceData(activeDomain_.mesh(), *core.faceShear)));
    }
    if (core.distances) {
      result.distances = std::array<std::vector<double>, 2>{
          vertexDataToVector(subdivision.interpolateAcrossB(
              activeVertexData(activeDomain_.mesh(), (*core.distances)[0]))),
          vertexDataToVector(subdivision.interpolateAcrossB(
              activeVertexData(activeDomain_.mesh(), (*core.distances)[1]))),
      };
    }
    if (sampleVertexShear && core.faceShear) {
      result.vertexShear =
          averageFaceScalarsToVertices(
              SurfaceMeshData{*result.mesh.vertices3D, result.mesh.faces},
              *result.faceShear,
              FaceScalarAveraging::FaceArea);
    }
  } else {
    result.mesh = makeExtrinsicResultMesh(reference_.meshData());
    result.origin.intrinsicPoint = core.intrinsicSeed;
    result.origin.intrinsicFamilyDirections = {core.intrinsicDirections[0], core.intrinsicDirections[2]};
    result.origin.extrinsicPoint = extrinsicSeed->cartesian;
    result.origin.extrinsicFamilyDirections = {core.cartesianDirections[0], core.cartesianDirections[2]};
    result.traces = makeTraceFamilies(core.generators, ResultDomain::Extrinsic);

    if (core.distances) {
      if ((*core.distances)[0].size() == reference_.meshData().vertices.size()) {
        result.distances = core.distances;
      } else {
        result.distances = std::array<std::vector<double>, 2>{
            restrictVertexScalarsToInput(activeDomain_.triangulation(), (*core.distances)[0]),
            restrictVertexScalarsToInput(activeDomain_.triangulation(), (*core.distances)[1]),
        };
      }
    }

    if (sampleVertexShear && core.faceShear) {
      const std::vector<double> activeVertexShear =
          averageIntrinsicFaceScalarsToVertices(activeDomain_.mesh(), *core.faceShear);
      if (activeVertexShear.size() == reference_.meshData().vertices.size()) {
        result.vertexShear = activeVertexShear;
      } else {
        result.vertexShear = restrictVertexScalarsToInput(activeDomain_.triangulation(), activeVertexShear);
      }
    }
  }

  return result;
}

} // namespace geodesic_draping

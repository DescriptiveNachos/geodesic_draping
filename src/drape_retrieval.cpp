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

DrapeTrace makeTrace(gcs::VertexPositionGeometry& inputGeometry,
                     gcs::IntrinsicTriangulation& triangulation,
                     const IntrinsicGeneratorTrace& generator,
                     ResultDomain domain) {
  DrapeTrace trace;
  trace.hitBoundary = generator.hitBoundary;
  trace.length = generator.length;
  if (domain == ResultDomain::Intrinsic) {
    trace.intrinsicPoints.reserve(generator.points.size());
    for (const gcs::SurfacePoint& point : generator.points) {
      trace.intrinsicPoints.push_back(toSurfaceReference(point));
    }
  } else {
    trace.extrinsicPoints.reserve(generator.points.size());
    for (const gcs::SurfacePoint& point : generator.points) {
      trace.extrinsicPoints.push_back(
          interpolateSurfacePoint(triangulation.equivalentPointOnInput(point), inputGeometry));
    }
  }
  return trace;
}

} // namespace

SeedProjection toExtrinsicSeed(gcs::VertexPositionGeometry& inputGeometry,
                               gcs::IntrinsicTriangulation& triangulation,
                               const gcs::SurfacePoint& intrinsicSeed) {
  const gcs::SurfacePoint inputPoint =
      triangulation.equivalentPointOnInput(intrinsicSeed).inSomeFace();

  SeedProjection seed;
  seed.surfacePoint.faceIndex = inputPoint.face.getIndex();
  seed.surfacePoint.barycentric =
      Vec3(inputPoint.faceCoords.x, inputPoint.faceCoords.y, inputPoint.faceCoords.z);
  seed.cartesian = interpolateSurfacePoint(inputPoint, inputGeometry);
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

std::array<TraceFamily, 2> makeTraceFamilies(
    gcs::VertexPositionGeometry& inputGeometry,
    gcs::IntrinsicTriangulation& triangulation,
    const std::array<IntrinsicGeneratorTrace, 4>& generators,
    ResultDomain domain) {
  return {
      TraceFamily{
          makeTrace(inputGeometry, triangulation, generators[0], domain),
          makeTrace(inputGeometry, triangulation, generators[1], domain)},
      TraceFamily{
          makeTrace(inputGeometry, triangulation, generators[2], domain),
          makeTrace(inputGeometry, triangulation, generators[3], domain)},
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
          ? std::optional<SeedProjection>{toExtrinsicSeed(
                *inputSurface_.geometry,
                *intrinsicTriangulation_,
                core.intrinsicSeed)}
          : std::nullopt;
  const SurfaceReference intrinsicSeed = toSurfaceReference(core.intrinsicSeed);
  const std::array<TangentVectorRef, 4> intrinsicDirections = {
      toTangentVectorRef(core.intrinsicDirections[0]),
      toTangentVectorRef(core.intrinsicDirections[1]),
      toTangentVectorRef(core.intrinsicDirections[2]),
      toTangentVectorRef(core.intrinsicDirections[3]),
  };
  const std::array<Vec3, 4> extrinsicDirections = cartesianDirectionsFromIntrinsic(
      *inputSurface_.geometry,
      *intrinsicTriangulation_,
      core.intrinsicSeed,
      core.intrinsicDirections);

  if (retrieval == ResultDomain::Intrinsic) {
    result.mesh = makeIntrinsicResultMesh(*intrinsicTriangulation_->intrinsicMesh, *intrinsicTriangulation_);
    result.origin.intrinsicPoint = intrinsicSeed;
    result.origin.intrinsicFamilyDirections = {intrinsicDirections[0], intrinsicDirections[2]};
    result.traces = makeTraceFamilies(
        *inputSurface_.geometry,
        *intrinsicTriangulation_,
        core.generators,
        ResultDomain::Intrinsic);
    result.directions = core.directions;
    result.faceShear = core.faceShear;
    result.distances = core.distances;
    if (sampleVertexShear && core.faceShear) {
      result.vertexShear =
          averageIntrinsicFaceScalarsToVertices(*intrinsicTriangulation_->intrinsicMesh, *core.faceShear);
    }
  } else if (retrieval == ResultDomain::Subdivision) {
    gcs::CommonSubdivision& subdivision = intrinsicTriangulation_->getCommonSubdivision();
    subdivision.constructMesh();
    result.mesh = makeSubdivisionResultMesh(subdivision, *inputSurface_.geometry);
    result.mesh.intrinsicSourceFaceColor =
        faceDataToVector(subdivision.copyFromB(gcs::niceColors(*intrinsicTriangulation_->intrinsicMesh)));
    result.origin.intrinsicPoint = intrinsicSeed;
    result.origin.intrinsicFamilyDirections = {intrinsicDirections[0], intrinsicDirections[2]};
    result.origin.extrinsicPoint = extrinsicSeed->cartesian;
    result.origin.extrinsicFamilyDirections = {extrinsicDirections[0], extrinsicDirections[2]};
    result.traces = makeTraceFamilies(
        *inputSurface_.geometry,
        *intrinsicTriangulation_,
        core.generators,
        ResultDomain::Subdivision);

    result.directions = {
        faceVectorDataToVector(subdivision.copyFromB(
            activeFaceVectorData(*intrinsicTriangulation_->intrinsicMesh, core.directions[0]))),
        faceVectorDataToVector(subdivision.copyFromB(
            activeFaceVectorData(*intrinsicTriangulation_->intrinsicMesh, core.directions[1]))),
    };
    if (core.faceShear) {
      result.faceShear = faceDataToVector(subdivision.copyFromB(
          activeFaceData(*intrinsicTriangulation_->intrinsicMesh, *core.faceShear)));
    }
    if (core.distances) {
      result.distances = std::array<std::vector<double>, 2>{
          vertexDataToVector(subdivision.interpolateAcrossB(
              activeVertexData(*intrinsicTriangulation_->intrinsicMesh, (*core.distances)[0]))),
          vertexDataToVector(subdivision.interpolateAcrossB(
              activeVertexData(*intrinsicTriangulation_->intrinsicMesh, (*core.distances)[1]))),
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
    result.mesh = makeExtrinsicResultMesh(meshData_);
    result.origin.intrinsicPoint = intrinsicSeed;
    result.origin.intrinsicFamilyDirections = {intrinsicDirections[0], intrinsicDirections[2]};
    result.origin.extrinsicPoint = extrinsicSeed->cartesian;
    result.origin.extrinsicFamilyDirections = {extrinsicDirections[0], extrinsicDirections[2]};
    result.traces = makeTraceFamilies(
        *inputSurface_.geometry,
        *intrinsicTriangulation_,
        core.generators,
        ResultDomain::Extrinsic);

    if (core.distances) {
      if ((*core.distances)[0].size() == meshData_.vertices.size()) {
        result.distances = core.distances;
      } else {
        result.distances = std::array<std::vector<double>, 2>{
            restrictVertexScalarsToInput(*intrinsicTriangulation_, (*core.distances)[0]),
            restrictVertexScalarsToInput(*intrinsicTriangulation_, (*core.distances)[1]),
        };
      }
    }

    if (sampleVertexShear && core.faceShear) {
      const std::vector<double> activeVertexShear =
          averageIntrinsicFaceScalarsToVertices(*intrinsicTriangulation_->intrinsicMesh, *core.faceShear);
      if (activeVertexShear.size() == meshData_.vertices.size()) {
        result.vertexShear = activeVertexShear;
      } else {
        result.vertexShear = restrictVertexScalarsToInput(*intrinsicTriangulation_, activeVertexShear);
      }
    }
  }

  return result;
}

} // namespace geodesic_draping

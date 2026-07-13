#include "geodrape_internal.h"

#include <limits>
#include <tuple>

namespace geodesic_draping {
namespace {

bool isVertexSurfacePoint(const gcs::SurfacePoint& point, gcs::Vertex vertex) {
  return point.type == gcs::SurfacePointType::Vertex && point.vertex == vertex;
}

std::optional<double> parameterAlongEdge(const gcs::SurfacePoint& point, gcs::Edge edge) {
  if (point.type == gcs::SurfacePointType::Vertex) {
    if (point.vertex == edge.halfedge().tailVertex()) {
      return 0.0;
    }
    if (point.vertex == edge.halfedge().tipVertex()) {
      return 1.0;
    }
    return std::nullopt;
  }
  if (point.type == gcs::SurfacePointType::Edge && point.edge == edge) {
    return point.tEdge;
  }
  return std::nullopt;
}

} // namespace

CommonSubdivisionDebugInfo GeoDrapeSolver::debugCommonSubdivision(bool attemptMeshConstruction) {
  gcs::CommonSubdivision& subdivision = intrinsicTriangulation_->getCommonSubdivision();

  CommonSubdivisionDebugInfo info;
  info.rawSubdivisionPointCount = subdivision.subdivisionPoints.size();

  std::vector<bool> representedInputVertices(inputSurface_.mesh->nVertices(), false);
  for (const gcs::CommonSubdivisionPoint& point : subdivision.subdivisionPoints) {
    switch (point.intersectionType) {
    case gcs::CSIntersectionType::VERTEX_VERTEX:
      ++info.vertexVertexCount;
      if (point.posA.type == gcs::SurfacePointType::Vertex &&
          point.posB.type == gcs::SurfacePointType::Vertex &&
          point.posA.vertex.getIndex() < representedInputVertices.size()) {
        representedInputVertices[point.posA.vertex.getIndex()] = true;
      }
      break;
    case gcs::CSIntersectionType::EDGE_TRANSVERSE:
      ++info.edgeTransverseCount;
      break;
    case gcs::CSIntersectionType::EDGE_PARALLEL:
      ++info.edgeParallelCount;
      break;
    case gcs::CSIntersectionType::FACE_VERTEX:
      ++info.faceVertexCount;
      break;
    case gcs::CSIntersectionType::EDGE_VERTEX:
      ++info.edgeVertexCount;
      break;
    }
  }
  for (bool represented : representedInputVertices) {
    if (!represented) {
      ++info.missingInputVertexCount;
    }
  }

  auto inspectEdgeList = [](const std::vector<gcs::CommonSubdivisionPoint*>& points,
                            gcs::Edge edge,
                            bool inspectA,
                            size_t& emptyCount,
                            size_t& invalidEndpointCount,
                            size_t& nonMonotoneCount) {
    if (points.empty()) {
      ++emptyCount;
      return;
    }

    const gcs::SurfacePoint& first = inspectA ? points.front()->posA : points.front()->posB;
    const gcs::SurfacePoint& last = inspectA ? points.back()->posA : points.back()->posB;
    if (!isVertexSurfacePoint(first, edge.halfedge().tailVertex()) ||
        !isVertexSurfacePoint(last, edge.halfedge().tipVertex())) {
      ++invalidEndpointCount;
    }

    double previous = -std::numeric_limits<double>::infinity();
    for (const gcs::CommonSubdivisionPoint* point : points) {
      const gcs::SurfacePoint& surfacePoint = inspectA ? point->posA : point->posB;
      const std::optional<double> t = parameterAlongEdge(surfacePoint, edge);
      if (!t) {
        continue;
      }
      if (*t + 1e-10 < previous) {
        ++nonMonotoneCount;
        return;
      }
      previous = *t;
    }
  };

  for (gcs::Edge edge : inputSurface_.mesh->edges()) {
    inspectEdgeList(
        subdivision.pointsAlongA[edge],
        edge,
        true,
        info.emptyPointsAlongACount,
        info.invalidPointsAlongAEndpointCount,
        info.nonMonotonePointsAlongACount);
  }
  for (gcs::Edge edge : intrinsicTriangulation_->intrinsicMesh->edges()) {
    inspectEdgeList(
        subdivision.pointsAlongB[edge],
        edge,
        false,
        info.emptyPointsAlongBCount,
        info.invalidPointsAlongBEndpointCount,
        info.nonMonotonePointsAlongBCount);
  }

  std::tie(
      info.expectedConstructedVertexCount,
      info.expectedConstructedEdgeCount,
      info.expectedConstructedFaceCount) = subdivision.elementCounts();

  if (attemptMeshConstruction) {
    info.attemptedMeshConstruction = true;
    try {
      subdivision.constructMesh();
      info.meshConstructed = true;
      info.constructedVertexCount = subdivision.mesh->nVertices();
      info.constructedFaceCount = subdivision.mesh->nFaces();
    } catch (const std::exception& e) {
      info.constructionError = e.what();
    }
  }

  return info;
}

} // namespace geodesic_draping

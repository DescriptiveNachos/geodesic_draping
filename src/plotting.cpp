#include "geodesic_draping/plotting.h"

#include <stdexcept>

#if GEODESIC_DRAPING_HAS_POLYSCOPE
#include "polyscope/options.h"
#include "polyscope/curve_network.h"
#include "polyscope/point_cloud.h"
#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"

#include <array>
#include <vector>
#endif

namespace geodesic_draping {

#if GEODESIC_DRAPING_HAS_POLYSCOPE
namespace {

std::vector<std::array<double, 3>> toPolyscopePoints(const std::vector<Vec3>& points) {
  std::vector<std::array<double, 3>> out;
  out.reserve(points.size());
  for (const Vec3& point : points) {
    out.push_back({point.x(), point.y(), point.z()});
  }
  return out;
}

std::vector<std::array<double, 3>> toPolyscopePoints(const SurfaceMeshData& mesh) {
  std::vector<std::array<double, 3>> out;
  out.reserve(mesh.vertices.size());
  for (const Vec3& vertex : mesh.vertices) {
    out.push_back({vertex.x(), vertex.y(), vertex.z()});
  }
  return out;
}

std::vector<std::array<size_t, 3>> toPolyscopeFaces(const SurfaceMeshData& mesh) {
  std::vector<std::array<size_t, 3>> out;
  out.reserve(mesh.faces.size());
  for (const Face& face : mesh.faces) {
    out.push_back(face);
  }
  return out;
}

std::vector<Vec3> pairedGeneratorPoints(const GeneratorTrace& negativeTrace,
                                        const GeneratorTrace& positiveTrace) {
  std::vector<Vec3> points;
  points.reserve(negativeTrace.points.size() + positiveTrace.points.size());
  points.insert(points.end(), negativeTrace.points.rbegin(), negativeTrace.points.rend());
  points.insert(points.end(), positiveTrace.points.begin(), positiveTrace.points.end());
  return points;
}

void ensurePolyscopeReady(bool clearExisting) {
  if (!polyscope::isInitialized()) {
    polyscope::init();
  }
  if (clearExisting) {
    polyscope::removeAllStructures();
  }
}

} // namespace
#endif

void plotSeedProjectionStep(const SurfaceMeshData& mesh,
                            const SeedProjection& projection,
                            const std::array<Vec3, 4>& directions,
                            const ProjectionPlotOptions& options) {
#if GEODESIC_DRAPING_HAS_POLYSCOPE
  ensurePolyscopeReady(options.clearExisting);

  polyscope::registerSurfaceMesh(options.name + " mesh", toPolyscopePoints(mesh), toPolyscopeFaces(mesh));

  const Vec3& origin = projection.cartesian;
  polyscope::registerPointCloud(options.name + " origin",
                                std::vector<std::array<double, 3>>{{origin.x(), origin.y(), origin.z()}});

  const std::array<std::string, 4> labels = {"dir +0", "dir -0", "dir +90", "dir -90"};
  for (size_t i = 0; i < directions.size(); ++i) {
    const Vec3 endpoint = origin + options.directionLength * directions[i];
    polyscope::registerCurveNetworkLine(
        options.name + " " + labels[i],
        std::vector<std::array<double, 3>>{
            {origin.x(), origin.y(), origin.z()},
            {endpoint.x(), endpoint.y(), endpoint.z()},
        });
  }

  if (options.show) {
    polyscope::show();
  }
#else
  (void)mesh;
  (void)projection;
  (void)directions;
  (void)options;
  throw std::runtime_error(
      "Polyscope plotting is disabled. Reconfigure with GEODESIC_DRAPING_ENABLE_POLYSCOPE=ON.");
#endif
}

void plotGeneratorTraces(const std::array<GeneratorTrace, 4>& traces,
                         const ProjectionPlotOptions& options) {
#if GEODESIC_DRAPING_HAS_POLYSCOPE
  ensurePolyscopeReady(options.clearExisting);

  for (size_t i = 0; i < traces.size(); ++i) {
    polyscope::registerCurveNetworkLine(options.name + " generator " + std::to_string(i),
                                        toPolyscopePoints(traces[i].points));
  }

  if (options.show) {
    polyscope::show();
  }
#else
  (void)traces;
  (void)options;
  throw std::runtime_error(
      "Polyscope plotting is disabled. Reconfigure with GEODESIC_DRAPING_ENABLE_POLYSCOPE=ON.");
#endif
}

void plotCompleteDrapeResult(const SurfaceMeshData& mesh,
                             const CompleteDrapeResult& result,
                             const ProjectionPlotOptions& options) {
#if GEODESIC_DRAPING_HAS_POLYSCOPE
  ensurePolyscopeReady(options.clearExisting);

  polyscope::SurfaceMesh* psMesh =
      polyscope::registerSurfaceMesh(options.name + " mesh", toPolyscopePoints(mesh), toPolyscopeFaces(mesh));

  psMesh->addVertexSignedDistanceQuantity("dist_0", result.distances[0]);
  psMesh->addVertexSignedDistanceQuantity("dist_1", result.distances[1]);
  psMesh->addVertexVectorQuantity("grad_0", toPolyscopePoints(result.gradients[0]), polyscope::VectorType::AMBIENT);
  psMesh->addVertexVectorQuantity("grad_1", toPolyscopePoints(result.gradients[1]), polyscope::VectorType::AMBIENT);
  psMesh->addVertexScalarQuantity("shear_degrees", result.shearAnglesDegrees)->setEnabled(true);

  const Vec3& origin = result.seed.cartesian;
  polyscope::registerPointCloud(options.name + " origin",
                                std::vector<std::array<double, 3>>{{origin.x(), origin.y(), origin.z()}});

  const std::array<std::string, 4> labels = {"dir +0", "dir -0", "dir +90", "dir -90"};
  for (size_t i = 0; i < result.directions.size(); ++i) {
    const Vec3 endpoint = origin + options.directionLength * result.directions[i];
    polyscope::registerCurveNetworkLine(
        options.name + " " + labels[i],
        std::vector<std::array<double, 3>>{
            {origin.x(), origin.y(), origin.z()},
            {endpoint.x(), endpoint.y(), endpoint.z()},
        });
  }

  for (size_t i = 0; i < result.generators.size(); ++i) {
    polyscope::registerCurveNetworkLine(options.name + " generator " + std::to_string(i),
                                        toPolyscopePoints(result.generators[i].points));
  }

  for (size_t i = 0; i < 2; ++i) {
    const std::vector<Vec3> points =
        pairedGeneratorPoints(result.generators[2 * i + 1], result.generators[2 * i]);
    polyscope::registerCurveNetworkLine(options.name + " source curve " + std::to_string(i),
                                        toPolyscopePoints(points));
  }

  if (options.show) {
    polyscope::show();
  }
#else
  (void)mesh;
  (void)result;
  (void)options;
  throw std::runtime_error(
      "Polyscope plotting is disabled. Reconfigure with GEODESIC_DRAPING_ENABLE_POLYSCOPE=ON.");
#endif
}

} // namespace geodesic_draping

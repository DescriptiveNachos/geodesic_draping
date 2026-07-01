#include "geodesic_draping/plotting.h"

#include <stdexcept>

#if GEODESIC_DRAPING_HAS_POLYSCOPE
#include "polyscope/curve_network.h"
#include "polyscope/point_cloud.h"
#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"

#include <array>
#include <vector>
#endif

namespace geodesic_draping {

void plotSeedProjectionStep(const SurfaceMeshData& mesh,
                            const SeedProjection& projection,
                            const std::array<Vec3, 4>& directions,
                            const ProjectionPlotOptions& options) {
#if GEODESIC_DRAPING_HAS_POLYSCOPE
  if (!polyscope::isInitialized()) {
    polyscope::init();
  }
  if (options.clearExisting) {
    polyscope::removeAllStructures();
  }

  std::vector<std::array<double, 3>> vertices;
  vertices.reserve(mesh.vertices.size());
  for (const Vec3& vertex : mesh.vertices) {
    vertices.push_back({vertex.x(), vertex.y(), vertex.z()});
  }

  std::vector<std::array<size_t, 3>> faces;
  faces.reserve(mesh.faces.size());
  for (const Face& face : mesh.faces) {
    faces.push_back(face);
  }

  polyscope::registerSurfaceMesh(options.name + " mesh", vertices, faces);

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
  if (!polyscope::isInitialized()) {
    polyscope::init();
  }
  if (options.clearExisting) {
    polyscope::removeAllStructures();
  }

  for (size_t i = 0; i < traces.size(); ++i) {
    std::vector<std::array<double, 3>> points;
    points.reserve(traces[i].points.size());
    for (const Vec3& point : traces[i].points) {
      points.push_back({point.x(), point.y(), point.z()});
    }
    polyscope::registerCurveNetworkLine(options.name + " generator " + std::to_string(i), points);
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

} // namespace geodesic_draping

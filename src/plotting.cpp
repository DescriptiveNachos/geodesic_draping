#include "geodesic_draping/plotting.h"

#include "geodesic_draping/diagnostics.h"

#include <stdexcept>

#if GEODESIC_DRAPING_HAS_POLYSCOPE
#include "polyscope/options.h"
#include "polyscope/curve_network.h"
#include "polyscope/point_cloud.h"
#include "polyscope/polyscope.h"
#include "polyscope/quantity.h"
#include "polyscope/surface_mesh.h"

#include <array>
#include <memory>
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
    polyscope::state::userCallback = nullptr;
  }
}

struct DrapeComparisonQuantities {
  polyscope::Quantity* completeDistance0 = nullptr;
  polyscope::Quantity* completeDistance1 = nullptr;
  polyscope::Quantity* completeGradient0 = nullptr;
  polyscope::Quantity* completeGradient1 = nullptr;
  polyscope::Quantity* completeGradient0Deviation = nullptr;
  polyscope::Quantity* completeGradient1Deviation = nullptr;
  polyscope::Quantity* completeShear = nullptr;
  polyscope::Quantity* fastGradient0 = nullptr;
  polyscope::Quantity* fastGradient1 = nullptr;
  polyscope::Quantity* fastGradient0Deviation = nullptr;
  polyscope::Quantity* fastGradient1Deviation = nullptr;
  polyscope::Quantity* fastShear = nullptr;
  bool showFast = false;
  bool showVectors = false;
  int completeScalar = 0;
};

void setEnabled(polyscope::Quantity* quantity, bool enabled) {
  if (quantity != nullptr) {
    quantity->setEnabled(enabled);
  }
}

void applyDrapeComparisonDisplay(DrapeComparisonQuantities& quantities) {
  setEnabled(quantities.completeDistance0, false);
  setEnabled(quantities.completeDistance1, false);
  setEnabled(quantities.completeGradient0, false);
  setEnabled(quantities.completeGradient1, false);
  setEnabled(quantities.completeGradient0Deviation, false);
  setEnabled(quantities.completeGradient1Deviation, false);
  setEnabled(quantities.completeShear, false);
  setEnabled(quantities.fastGradient0, false);
  setEnabled(quantities.fastGradient1, false);
  setEnabled(quantities.fastGradient0Deviation, false);
  setEnabled(quantities.fastGradient1Deviation, false);
  setEnabled(quantities.fastShear, false);

  if (quantities.showFast) {
    setEnabled(quantities.fastShear, true);
    setEnabled(quantities.fastGradient0, quantities.showVectors);
    setEnabled(quantities.fastGradient1, quantities.showVectors);
    setEnabled(quantities.fastGradient0Deviation, true);
    setEnabled(quantities.fastGradient1Deviation, true);
  } else {
    setEnabled(quantities.completeShear, quantities.completeScalar == 0);
    setEnabled(quantities.completeDistance0, quantities.completeScalar == 1);
    setEnabled(quantities.completeDistance1, quantities.completeScalar == 2);
    setEnabled(quantities.completeGradient0, quantities.showVectors);
    setEnabled(quantities.completeGradient1, quantities.showVectors);
    setEnabled(quantities.completeGradient0Deviation, true);
    setEnabled(quantities.completeGradient1Deviation, true);
  }
}

void registerSeedDirectionsAndGenerators(const std::string& name,
                                         const SeedProjection& seed,
                                         const std::array<Vec3, 4>& directions,
                                         const std::array<GeneratorTrace, 4>& generators,
                                         double directionLength) {
  const Vec3& origin = seed.cartesian;
  polyscope::registerPointCloud(name + " origin",
                                std::vector<std::array<double, 3>>{{origin.x(), origin.y(), origin.z()}});

  const std::array<std::string, 4> labels = {"dir +0", "dir -0", "dir +90", "dir -90"};
  for (size_t i = 0; i < directions.size(); ++i) {
    const Vec3 endpoint = origin + directionLength * directions[i];
    polyscope::registerCurveNetworkLine(
        name + " " + labels[i],
        std::vector<std::array<double, 3>>{
            {origin.x(), origin.y(), origin.z()},
            {endpoint.x(), endpoint.y(), endpoint.z()},
        });
  }

  for (size_t i = 0; i < generators.size(); ++i) {
    polyscope::registerCurveNetworkLine(name + " generator " + std::to_string(i),
                                        toPolyscopePoints(generators[i].points));
  }

  for (size_t i = 0; i < 2; ++i) {
    const std::vector<Vec3> points =
        pairedGeneratorPoints(generators[2 * i + 1], generators[2 * i]);
    polyscope::registerCurveNetworkLine(name + " source curve " + std::to_string(i),
                                        toPolyscopePoints(points));
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

  registerSeedDirectionsAndGenerators(
      options.name, result.seed, result.directions, result.generators, options.directionLength);

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

void plotDrapeComparisonResult(const SurfaceMeshData& mesh,
                               const CompleteDrapeResult& completeResult,
                               const FastDrapeResult& fastResult,
                               const ProjectionPlotOptions& options) {
#if GEODESIC_DRAPING_HAS_POLYSCOPE
  ensurePolyscopeReady(options.clearExisting);

  polyscope::SurfaceMesh* psMesh =
      polyscope::registerSurfaceMesh(options.name + " mesh", toPolyscopePoints(mesh), toPolyscopeFaces(mesh));

  auto quantities = std::make_shared<DrapeComparisonQuantities>();
  quantities->completeDistance0 = psMesh->addVertexSignedDistanceQuantity(
      "complete dist_0", completeResult.distances[0]);
  quantities->completeDistance1 = psMesh->addVertexSignedDistanceQuantity(
      "complete dist_1", completeResult.distances[1]);
  quantities->completeGradient0 = psMesh->addVertexVectorQuantity(
      "complete grad_0", toPolyscopePoints(completeResult.gradients[0]), polyscope::VectorType::AMBIENT);
  quantities->completeGradient1 = psMesh->addVertexVectorQuantity(
      "complete grad_1", toPolyscopePoints(completeResult.gradients[1]), polyscope::VectorType::AMBIENT);
  const auto completeMagnitudeDiagnostics = analyzeCompleteGradientMagnitudes(completeResult);
  quantities->completeGradient0Deviation = psMesh->addVertexScalarQuantity(
      "abs(|complete grad_0| - 1)", completeMagnitudeDiagnostics[0].absDeviationFromUnit);
  quantities->completeGradient1Deviation = psMesh->addVertexScalarQuantity(
      "abs(|complete grad_1| - 1)", completeMagnitudeDiagnostics[1].absDeviationFromUnit);
  quantities->completeShear = psMesh->addVertexScalarQuantity(
      "complete shear_degrees", completeResult.shearAnglesDegrees);
  quantities->fastGradient0 = psMesh->addVertexVectorQuantity(
      "fast grad_0", toPolyscopePoints(fastResult.gradients[0]), polyscope::VectorType::AMBIENT);
  quantities->fastGradient1 = psMesh->addVertexVectorQuantity(
      "fast grad_1", toPolyscopePoints(fastResult.gradients[1]), polyscope::VectorType::AMBIENT);
  const auto fastMagnitudeDiagnostics = analyzeFastGradientMagnitudes(fastResult);
  quantities->fastGradient0Deviation = psMesh->addVertexScalarQuantity(
      "abs(|fast grad_0| - 1)", fastMagnitudeDiagnostics[0].absDeviationFromUnit);
  quantities->fastGradient1Deviation = psMesh->addVertexScalarQuantity(
      "abs(|fast grad_1| - 1)", fastMagnitudeDiagnostics[1].absDeviationFromUnit);
  quantities->fastShear = psMesh->addVertexScalarQuantity(
      "fast shear_degrees", fastResult.shearAnglesDegrees);
  applyDrapeComparisonDisplay(*quantities);

  registerSeedDirectionsAndGenerators(options.name,
                                      completeResult.seed,
                                      completeResult.directions,
                                      completeResult.generators,
                                      options.directionLength);

  polyscope::state::userCallback = [quantities]() {
    bool changed = false;
    changed |= ImGui::Checkbox("Show fast result", &quantities->showFast);
    changed |= ImGui::Checkbox("Show vector fields", &quantities->showVectors);

    if (!quantities->showFast) {
      ImGui::TextUnformatted("Complete scalar field");
      changed |= ImGui::RadioButton("shear", &quantities->completeScalar, 0);
      changed |= ImGui::RadioButton("dist_0", &quantities->completeScalar, 1);
      changed |= ImGui::RadioButton("dist_1", &quantities->completeScalar, 2);
    } else {
      ImGui::TextUnformatted("Fast mode displays fast shear");
    }

    if (changed) {
      applyDrapeComparisonDisplay(*quantities);
      polyscope::requestRedraw();
    }
  };

  if (options.show) {
    polyscope::show();
  }
#else
  (void)mesh;
  (void)completeResult;
  (void)fastResult;
  (void)options;
  throw std::runtime_error(
      "Polyscope plotting is disabled. Reconfigure with GEODESIC_DRAPING_ENABLE_POLYSCOPE=ON.");
#endif
}

} // namespace geodesic_draping

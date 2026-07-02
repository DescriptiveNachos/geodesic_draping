#include "geodesic_draping/plotting.h"

#include "geodesic_draping/diagnostics.h"
#include "geodesic_draping/field_processing.h"
#include "geodesic_draping/geometrycentral_adapter.h"

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
  polyscope::Quantity* fastFaceDirection0 = nullptr;
  polyscope::Quantity* fastFaceDirection1 = nullptr;
  polyscope::Quantity* fastFaceDirection0Deviation = nullptr;
  polyscope::Quantity* fastFaceDirection1Deviation = nullptr;
  polyscope::Quantity* fastFaceShear = nullptr;
  polyscope::Quantity* fastVertexShear = nullptr;
  bool showFast = false;
  bool showVectors = false;
  bool showFaceVectors = false;
  int completeScalar = 0;
  int fastScalar = 0;
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
  setEnabled(quantities.fastFaceDirection0, false);
  setEnabled(quantities.fastFaceDirection1, false);
  setEnabled(quantities.fastFaceDirection0Deviation, false);
  setEnabled(quantities.fastFaceDirection1Deviation, false);
  setEnabled(quantities.fastFaceShear, false);
  setEnabled(quantities.fastVertexShear, false);

  if (quantities.showFast) {
    setEnabled(quantities.fastFaceShear, quantities.fastScalar == 0);
    setEnabled(quantities.fastVertexShear, quantities.fastScalar == 1);
    setEnabled(quantities.fastFaceDirection0, quantities.showFaceVectors);
    setEnabled(quantities.fastFaceDirection1, quantities.showFaceVectors);
    setEnabled(quantities.fastFaceDirection0Deviation, true);
    setEnabled(quantities.fastFaceDirection1Deviation, true);
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

void plotDrapeComparisonResult(const SurfaceMeshData& mesh,
                               const DrapeResult& completeResult,
                               const DrapeResult& fastResult,
                               const ProjectionPlotOptions& options) {
#if GEODESIC_DRAPING_HAS_POLYSCOPE
  if (!completeResult.distances || !completeResult.gradients || !completeResult.vertexShearAnglesDegrees) {
    throw std::runtime_error("plotDrapeComparisonResult requires a complete-mode result");
  }
  if (!fastResult.faceShearAnglesDegrees || !fastResult.vertexShearAnglesDegrees) {
    throw std::runtime_error("plotDrapeComparisonResult requires a fast result with sampled secondary shear");
  }

  ensurePolyscopeReady(options.clearExisting);

  polyscope::SurfaceMesh* psMesh =
      polyscope::registerSurfaceMesh(options.name + " mesh", toPolyscopePoints(mesh), toPolyscopeFaces(mesh));

  auto quantities = std::make_shared<DrapeComparisonQuantities>();
  quantities->completeDistance0 = psMesh->addVertexSignedDistanceQuantity(
      "complete dist_0", (*completeResult.distances)[0]);
  quantities->completeDistance1 = psMesh->addVertexSignedDistanceQuantity(
      "complete dist_1", (*completeResult.distances)[1]);
  quantities->completeGradient0 = psMesh->addVertexVectorQuantity(
      "complete grad_0", toPolyscopePoints((*completeResult.gradients)[0]), polyscope::VectorType::AMBIENT);
  quantities->completeGradient1 = psMesh->addVertexVectorQuantity(
      "complete grad_1", toPolyscopePoints((*completeResult.gradients)[1]), polyscope::VectorType::AMBIENT);
  const auto completeMagnitudeDiagnostics = analyzeGradientMagnitudes(*completeResult.gradients);
  quantities->completeGradient0Deviation = psMesh->addVertexScalarQuantity(
      "abs(|complete grad_0| - 1)", completeMagnitudeDiagnostics[0].absDeviationFromUnit);
  quantities->completeGradient1Deviation = psMesh->addVertexScalarQuantity(
      "abs(|complete grad_1| - 1)", completeMagnitudeDiagnostics[1].absDeviationFromUnit);
  quantities->completeShear = psMesh->addVertexScalarQuantity(
      "complete shear_degrees", *completeResult.vertexShearAnglesDegrees);
  GeometryCentralSurface surface = makeGeometryCentralSurface(mesh);
  const std::vector<Vec3> fastFaceDirections0 = faceDirectionsToExtrinsicVectors(
      surface,
      fastResult.faceDirections[0]);
  const std::vector<Vec3> fastFaceDirections1 = faceDirectionsToExtrinsicVectors(
      surface,
      fastResult.faceDirections[1]);
  const auto fastFaceMagnitudeDiagnostics0 = analyzeVectorMagnitudes(fastFaceDirections0);
  const auto fastFaceMagnitudeDiagnostics1 = analyzeVectorMagnitudes(fastFaceDirections1);
  quantities->fastFaceDirection0 = psMesh->addFaceVectorQuantity(
      "fast face dir_0", toPolyscopePoints(fastFaceDirections0), polyscope::VectorType::AMBIENT);
  quantities->fastFaceDirection1 = psMesh->addFaceVectorQuantity(
      "fast face dir_1", toPolyscopePoints(fastFaceDirections1), polyscope::VectorType::AMBIENT);
  quantities->fastFaceDirection0Deviation = psMesh->addFaceScalarQuantity(
      "abs(|fast face dir_0| - 1)", fastFaceMagnitudeDiagnostics0.absDeviationFromUnit);
  quantities->fastFaceDirection1Deviation = psMesh->addFaceScalarQuantity(
      "abs(|fast face dir_1| - 1)", fastFaceMagnitudeDiagnostics1.absDeviationFromUnit);
  quantities->fastFaceShear = psMesh->addFaceScalarQuantity(
      "fast face shear_degrees", *fastResult.faceShearAnglesDegrees);
  quantities->fastVertexShear = psMesh->addVertexScalarQuantity(
      "fast vertex shear_degrees", *fastResult.vertexShearAnglesDegrees);
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
    if (quantities->showFast) {
      changed |= ImGui::Checkbox("Show face vectors", &quantities->showFaceVectors);
    }

    if (!quantities->showFast) {
      ImGui::TextUnformatted("Complete scalar field");
      changed |= ImGui::RadioButton("shear", &quantities->completeScalar, 0);
      changed |= ImGui::RadioButton("dist_0", &quantities->completeScalar, 1);
      changed |= ImGui::RadioButton("dist_1", &quantities->completeScalar, 2);
    } else {
      ImGui::TextUnformatted("Fast scalar field");
      changed |= ImGui::RadioButton("face shear", &quantities->fastScalar, 0);
      changed |= ImGui::RadioButton("vertex shear", &quantities->fastScalar, 1);
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

#pragma once

#include "geodesic_draping/geodrape.h"
#include "geodesic_draping/mesh.h"

#include <array>
#include <cstddef>
#include <vector>

namespace geodesic_draping {

struct MagnitudeStats {
  double min = 0.0;
  double max = 0.0;
  double mean = 0.0;
  double maxAbsDeviationFromUnit = 0.0;
  double meanAbsDeviationFromUnit = 0.0;
  size_t finiteCount = 0;
  size_t nonFiniteCount = 0;
  size_t nearZeroCount = 0;
};

struct VectorMagnitudeDiagnostics {
  std::vector<double> magnitudes;
  std::vector<double> absDeviationFromUnit;
  MagnitudeStats stats;
};

VectorMagnitudeDiagnostics analyzeVectorMagnitudes(const std::vector<Vec3>& vectors,
                                                   double nearZeroEpsilon = 1e-12);

std::array<VectorMagnitudeDiagnostics, 2> analyzeCompleteGradientMagnitudes(
    const CompleteDrapeResult& result,
    double nearZeroEpsilon = 1e-12);

} // namespace geodesic_draping

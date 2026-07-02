#include "geodesic_draping/diagnostics.h"

#include <algorithm>
#include <cmath>
#include <limits>

namespace geodesic_draping {

VectorMagnitudeDiagnostics analyzeVectorMagnitudes(const std::vector<Vec3>& vectors,
                                                   double nearZeroEpsilon) {
  VectorMagnitudeDiagnostics diagnostics;
  diagnostics.magnitudes.reserve(vectors.size());
  diagnostics.absDeviationFromUnit.reserve(vectors.size());

  double magnitudeSum = 0.0;
  double deviationSum = 0.0;
  diagnostics.stats.min = std::numeric_limits<double>::infinity();
  diagnostics.stats.max = -std::numeric_limits<double>::infinity();

  for (const Vec3& vector : vectors) {
    const double magnitude = vector.norm();
    diagnostics.magnitudes.push_back(magnitude);

    if (!std::isfinite(magnitude)) {
      diagnostics.absDeviationFromUnit.push_back(std::numeric_limits<double>::quiet_NaN());
      ++diagnostics.stats.nonFiniteCount;
      continue;
    }

    const double deviation = std::abs(magnitude - 1.0);
    diagnostics.absDeviationFromUnit.push_back(deviation);
    diagnostics.stats.min = std::min(diagnostics.stats.min, magnitude);
    diagnostics.stats.max = std::max(diagnostics.stats.max, magnitude);
    diagnostics.stats.maxAbsDeviationFromUnit =
        std::max(diagnostics.stats.maxAbsDeviationFromUnit, deviation);
    magnitudeSum += magnitude;
    deviationSum += deviation;
    ++diagnostics.stats.finiteCount;
    if (magnitude <= nearZeroEpsilon) {
      ++diagnostics.stats.nearZeroCount;
    }
  }

  if (diagnostics.stats.finiteCount == 0) {
    diagnostics.stats.min = std::numeric_limits<double>::quiet_NaN();
    diagnostics.stats.max = std::numeric_limits<double>::quiet_NaN();
    diagnostics.stats.mean = std::numeric_limits<double>::quiet_NaN();
    diagnostics.stats.meanAbsDeviationFromUnit = std::numeric_limits<double>::quiet_NaN();
    return diagnostics;
  }

  const double finiteCount = static_cast<double>(diagnostics.stats.finiteCount);
  diagnostics.stats.mean = magnitudeSum / finiteCount;
  diagnostics.stats.meanAbsDeviationFromUnit = deviationSum / finiteCount;
  return diagnostics;
}

std::array<VectorMagnitudeDiagnostics, 2> analyzeCompleteGradientMagnitudes(
    const CompleteDrapeResult& result,
    double nearZeroEpsilon) {
  return {
      analyzeVectorMagnitudes(result.gradients[0], nearZeroEpsilon),
      analyzeVectorMagnitudes(result.gradients[1], nearZeroEpsilon),
  };
}

} // namespace geodesic_draping

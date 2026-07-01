#include "geodesic_draping/field_processing.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace geodesic_draping {

std::vector<double> computeShearAnglesDegrees(const std::vector<Vec3>& gradients0,
                                              const std::vector<Vec3>& gradients1) {
  if (gradients0.size() != gradients1.size()) {
    throw std::runtime_error("computeShearAnglesDegrees requires equally sized gradient arrays");
  }

  constexpr double pi = 3.141592653589793238462643383279502884;
  std::vector<double> shear;
  shear.reserve(gradients0.size());
  for (size_t i = 0; i < gradients0.size(); ++i) {
    const double norm0 = gradients0[i].norm();
    const double norm1 = gradients1[i].norm();
    const double cosTheta = std::clamp(gradients1[i].dot(gradients0[i]) / (norm0 * norm1), -1.0, 1.0);
    shear.push_back(std::abs((std::acos(cosTheta) * 180.0 / pi) - 90.0));
  }
  return shear;
}

} // namespace geodesic_draping

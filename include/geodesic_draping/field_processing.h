#pragma once

#include "geodesic_draping/mesh.h"

#include <vector>

namespace geodesic_draping {

std::vector<double> computeShearAnglesDegrees(const std::vector<Vec3>& gradients0,
                                              const std::vector<Vec3>& gradients1);

} // namespace geodesic_draping

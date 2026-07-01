#include <Eigen/Dense>

#include <cassert>
#include <cmath>

namespace {

bool near(double a, double b, double tolerance = 1e-12) {
  return std::abs(a - b) <= tolerance;
}

} // namespace

int main() {
  Eigen::Vector3d x(1.0, 0.0, 0.0);
  Eigen::Vector3d y(0.0, 1.0, 0.0);
  Eigen::Vector3d z = x.cross(y);

  assert(near(x.dot(y), 0.0));
  assert(near(z.x(), 0.0));
  assert(near(z.y(), 0.0));
  assert(near(z.z(), 1.0));
  assert(near(z.norm(), 1.0));

  Eigen::Matrix3d identity = Eigen::Matrix3d::Identity();
  Eigen::Vector3d transformed = identity * (x + y + z);
  assert(near(transformed.x(), 1.0));
  assert(near(transformed.y(), 1.0));
  assert(near(transformed.z(), 1.0));

  return 0;
}

#include "coordinate.h"

#include <math.h>

#include <random>

double Coordinate::length() const { return sqrt(x * x + y * y + z * z); }

double Coordinate::phi() const {
  double phi = atan2(y, x);
  if (phi < 0) {
    phi = 2.0 * M_PI + phi;
  }
  return phi;
}

Coordinate Coordinate::cross_product(const Coordinate& other) const {
  return Coordinate(y * other.z - z * other.y, z * other.x - x * other.z,
                    x * other.y - y * other.x);
}

double Coordinate::dot_product(const Coordinate& other) const {
  return x * other.x + y * other.y + z * other.z;
}

Coordinate Coordinate::coordinate_by_rotation(Coordinate& axis,
                                              double angle) const {
  /**
   * Rodrigues' formula consists of three summands.
   */

  // 1. Summand
  Coordinate first = *this * cos(angle);

  // 2. Summand
  Coordinate second = axis.cross_product(*this) * sin(angle);

  // 3. Summand
  Coordinate third = axis * axis.dot_product(*this) * (1.0 - cos(angle));

  return first + second + third;
}

Coordinate Coordinate::cartesian_to_spherical() const {
  double r = length();
  double lambda = atan2(y, x);
  if (lambda < 0) {
    lambda += 2.0 * M_PI;
  }
  double phi = atan2(z, sqrt(x * x + y * y));

  return Coordinate(r, lambda, phi);
}

Coordinate Coordinate::spherical_to_cartesian() const {
  return Coordinate(x * cos(z) * cos(y), x * cos(z) * sin(y), x * sin(z));
}

double Coordinate::central_angle_to_coordinate(const Coordinate& other) const {
  double lambda1 = this->y;
  double lambda2 = other.y;
  double phi1 = this->z;
  double phi2 = other.z;

  double deltaLambda = abs(lambda1 - lambda2);

  /**
   * The central angle is determined using the following function, which is said
   * to be rather robust against numerical difficulties.
   */
  return atan2(
      sqrt((cos(phi2) * sin(deltaLambda)) * (cos(phi2) * sin(deltaLambda)) +
           (cos(phi1) * sin(phi2) - sin(phi1) * cos(phi2) * cos(deltaLambda)) *
               (cos(phi1) * sin(phi2) -
                sin(phi1) * cos(phi2) * cos(deltaLambda))),
      (sin(phi1) * sin(phi2) + cos(phi1) * cos(phi2) * cos(deltaLambda)));
}

Coordinate Coordinate::normal_vector_for_plane_from_points(
    const Coordinate& coord1, const Coordinate& coord2,
    const Coordinate& coord3) {
  /**
   * We construct two vectors that lie in the plane.
   */
  if (coord1 == coord2 || coord1 == coord3 || coord2 == coord3) {
    return Coordinate();
  }

  Coordinate vector1 = coord2 - coord1;
  Coordinate vector2 = coord3 - coord1;

  /**
   * Now the cross product of these two vectors is orthogonal to them and is
   * therefore our normal vector.
   */

  Coordinate normal_vector = vector1.cross_product(vector2);

  /**
   * Due to numerical issues this can be zero even though the
   * coordinates are not really equal. Thus, we have to add this
   * check.
   */
  if (normal_vector.length() == 0.0) {
    return normal_vector;
  }

  return normal_vector / normal_vector.length();
}

Coordinate Coordinate::random_point_on_surface_of_unit_ball(int seed) {
  std::normal_distribution<> distribution(0.0, 1.0);
  std::default_random_engine generator(seed);
  double randomX = distribution(generator);
  double randomY = distribution(generator);
  double randomZ = distribution(generator);

  Coordinate candidate = Coordinate(randomX, randomY, randomZ);

  candidate /= candidate.length();

  return candidate;
}

/**
 * Operators
 */

Coordinate Coordinate::operator+(const Coordinate& other) const {
  return Coordinate(x + other.x, y + other.y, z + other.z);
}

Coordinate Coordinate::operator+=(const Coordinate& other) {
  x += other.x;
  y += other.y;
  z += other.z;
  return *this;
}

Coordinate Coordinate::operator-(const Coordinate& other) const {
  return Coordinate(x - other.x, y - other.y, z - other.z);
}

Coordinate Coordinate::operator*(const double factor) const {
  return Coordinate(factor * x, factor * y, factor * z);
}

Coordinate Coordinate::operator*=(const double factor) {
  x *= factor;
  y *= factor;
  z *= factor;
  return *this;
}

Coordinate Coordinate::operator/(const double factor) const {
  return Coordinate(x / factor, y / factor, z / factor);
}

Coordinate Coordinate::operator/=(const double factor) {
  x /= factor;
  y /= factor;
  z /= factor;
  return *this;
}

bool Coordinate::operator==(const Coordinate& other) const {
  return x == other.x && y == other.y && z == other.z;
}

#include "force.h"

#include <math.h>

#include "assert.h"

Force::Force() {}

Coordinate Force::force_from_coordinate_to_coordinate(
    const Coordinate &sender, const Coordinate &receiver, bool attractive) {
  assert(
      "The Force class is abstract. Derive from it and implement the "
      "'applyForceFromCoordinateToCoordinate' method!");

  return Coordinate();
}

double Force::radial_force_from_coordinate_to_coordinate(
    const Coordinate &sender, const Coordinate &receiver, bool attractive) {
  assert(
      "The Force class is abstract. Derive from it and implement the "
      "'radialForceFromCoordinateToCoordinate' method!");
  return 0.0;
}

double Force::distance_between_coordinates(const Coordinate &coord1,
                                           const Coordinate &coord2) const {
  assert(
      "The Force class is abstract. Derive from it and implement the "
      "'distanceBetweenCordinates' method!");
  return 0.0;
}

double Force::hyperbolic_distance_between_coordinates(
    const Coordinate &coord1, const Coordinate &coord2) {
  /**
   * Check whether one of the coordinates is the origin. If it is, we can simply
   * return the radius of the second coordinate.
   */
  Coordinate origin = Coordinate();
  if (coord1 == origin) {
    return coord2.length();
  }
  if (coord2 == origin) {
    return coord1.length();
  }

  /**
   * At first we determine central angle between the two vertices.
   * This is done by first converting the coordinates spherical coordinates.
   */

  Coordinate spherical1 = coord1.cartesian_to_spherical();
  Coordinate spherical2 = coord2.cartesian_to_spherical();

  double lambda1 = spherical1.y;
  double lambda2 = spherical2.y;
  double phi1 = spherical1.z;
  double phi2 = spherical2.z;

  if (lambda1 == lambda2 && phi1 == phi2) {
    return fabs(spherical1.x - spherical2.x);
  }

  double delta_lambda = fabs(lambda1 - lambda2);

  /**
   * The central angle is determined using the following function, which is said
   * to be rather robust against numerical difficulties.
   */
  double delta_sigma = atan2(
      sqrt((cos(phi2) * sin(delta_lambda)) * (cos(phi2) * sin(delta_lambda)) +
           (cos(phi1) * sin(phi2) - sin(phi1) * cos(phi2) * cos(delta_lambda)) *
               (cos(phi1) * sin(phi2) -
                sin(phi1) * cos(phi2) * cos(delta_lambda))),
      (sin(phi1) * sin(phi2) + cos(phi1) * cos(phi2) * cos(delta_lambda)));

  /**
   * Now we now the angular distance between the to particles. Next up: the
   * radial coordinates.
   */

  double r1 = spherical1.x;
  double r2 = spherical2.x;

  double distance =
      acosh(cosh(r1) * cosh(r2) - sinh(r1) * sinh(r2) * cos(delta_sigma));

  /**
   * We assume that NaNs occur when two coordinates are too close to
   * each other. In that case we set the distance to 0 and hope for
   * the best.
   */
  if (distance != distance) {
    distance = 0.0;
  }

  return distance;
}

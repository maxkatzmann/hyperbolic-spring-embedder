#include "constantforce.h"

#include <math.h>

ConstantForce::ConstantForce(double disk_radius, double maximal_force)
    : disk_radius(disk_radius), maximal_force(maximal_force){};

Coordinate ConstantForce::force_from_coordinate_to_coordinate(
    const Coordinate &sender, const Coordinate &receiver, bool attractive) {
  /**
   * The force depends on the distance between the two particles.
   */
  // double distance = distanceBetweenCoordinates(sender, receiver);
  /**
   * Depending on whether the force is attractive or not we determine the
   * coordinate of the receiver after the force is applied.
   */

  double rotation_angle = 0.0;

  /**
   * Now we determine the rotation axis. It is the normal vector of the
   * plane defined by the origin, the sender and the receiver.
   */
  Coordinate origin = Coordinate();
  Coordinate rotation_axis =
      Coordinate::normal_vector_for_plane_from_points(origin, receiver, sender);

  rotation_angle = 0.2 * maximal_force;

  if (!attractive) {
    rotation_axis *= -1;
  }
  /**
   * Applying the temperature.
   */
  rotation_angle *= temperature;

  return rotation_axis * rotation_angle;
}

double ConstantForce::distance_between_coordinates(
    const Coordinate &coord1, const Coordinate &coord2) const {
  return Force::hyperbolic_distance_between_coordinates(coord1, coord2);
}

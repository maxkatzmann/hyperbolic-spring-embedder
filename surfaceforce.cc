#include "surfaceforce.h"

#include <math.h>

SurfaceForce::SurfaceForce(double disk_radius, double maximal_force)
    : disk_radius(disk_radius), maximal_force(maximal_force){};

Coordinate SurfaceForce::force_from_coordinate_to_coordinate(
    const Coordinate &sender, const Coordinate &receiver, bool attractive) {
  /**
   * The force on the sphere surface is determined is simply a
   * rotation about an angle (the strength of the force) and a rotation
   * axis (the direction of the force)
   */
  double rotation_angle = 0.0;

  /**
   * First we determine the rotation axis. It is the normal vector of the
   * plane defined by the origin, the sender and the receiver.
   */
  Coordinate origin = Coordinate();
  Coordinate rotation_axis =
      Coordinate::normal_vector_for_plane_from_points(origin, receiver, sender);

  /**
   * If we're repulsive, we invert the axis of rotation.
   */
  if (!attractive) {
    rotation_axis *= -1.0;
  }

  /**
   * Determine the hyperbolic distance.
   */
  double hyperbolic_distance =
      Force::hyperbolic_distance_between_coordinates(sender, receiver);

  if (hyperbolic_distance == 0.0) {
    return Coordinate();
  }

  double hyperbolic_force = 0.0;

  if (attractive) {
    /**
     * Compute hyperbolic force.
     */
    double min_hyperbolic_distance = 0.5 * disk_radius;
    double max_hyperbolic_distance = 2.0 * disk_radius;

    if (hyperbolic_distance > min_hyperbolic_distance) {
      double hyperbolic_base =
          (hyperbolic_distance - min_hyperbolic_distance) /
          (max_hyperbolic_distance - min_hyperbolic_distance);

      hyperbolic_force = pow(hyperbolic_base, this->attractive_exponent);
    }
  } else {
    /**
     * Compute hyperbolic force.
     */
    double min_hyperbolic_distance = 1.0 * disk_radius;
    double max_hyperbolic_distance = 2.0 * disk_radius;

    if (hyperbolic_distance < min_hyperbolic_distance) {
      hyperbolic_force = 1.0;
    } else {
      double hyperbolic_base =
          (max_hyperbolic_distance - hyperbolic_distance) /
          (max_hyperbolic_distance - min_hyperbolic_distance);

      hyperbolic_force = pow(hyperbolic_base, this->repulsive_exponent);
    }
  }

  rotation_angle = hyperbolic_force;

  if (!attractive) {
    rotation_angle *= repulsive_force_coefficient;
  }

  /**
   * Until now the force values lie between 0 and 1. Now we scale them
   * with the maximal force and apply the temperature.
   */
  rotation_angle *= maximal_force * temperature;

  return rotation_axis * rotation_angle;
}

double SurfaceForce::radial_force_from_coordinate_to_coordinate(
    const Coordinate &sender, const Coordinate &receiver, bool attractive) {
  return 0.0;
}

double SurfaceForce::distance_between_coordinates(
    const Coordinate &coord1, const Coordinate &coord2) const {
  return Force::hyperbolic_distance_between_coordinates(coord1, coord2);
}

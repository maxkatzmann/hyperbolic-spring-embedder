//
//  surface.hpp
//  threege
//
//  Created by Maximilian Katzmann on 13.11.17.
//

#pragma once

#include <iostream>
#include <string>

#include "coordinate.h"
#include "force.h"

using namespace std;

class SurfaceForce : public Force {
 public:
  /**
   * The tangent force depends on the radius of the disk that the particles are
   * placed in.
   */
  double disk_radius;

  /**
   * Determines the maximum force that should be applied. I.e. the force that
   * is applied between the two vertices with the maximal distance.
   */
  double maximal_force;

  double repulsive_force_coefficient = 0.03;

  double attractive_exponent = 1.0;
  double repulsive_exponent = 0.5;

  SurfaceForce(double disk_radius, double maximal_force);

  SurfaceForce();

  Coordinate force_from_coordinate_to_coordinate(const Coordinate& sender,
                                                 const Coordinate& receiver,
                                                 bool attractive);

  double radial_force_from_coordinate_to_coordinate(const Coordinate& sender,
                                                    const Coordinate& receiver,
                                                    bool attractive);
  /**
   * The surface force uses an interpolation of the hyperbolic distance in the
   * 3-dimensional hyperbolic space and the distance on the sphere surface.
   *
   * @param  coord1 A coordinate.
   * @param  coord2 Another coordinate.
   * @return        The distance between the coordinates according to the
   * surface force.
   */
  double distance_between_coordinates(const Coordinate& coord1,
                                      const Coordinate& coord2) const;

  virtual ~SurfaceForce() = default;
  SurfaceForce(const SurfaceForce& other) = default;
  SurfaceForce(SurfaceForce&& other) = default;
  SurfaceForce& operator=(const SurfaceForce& other) = default;
  SurfaceForce& operator=(SurfaceForce&& other) = default;
};

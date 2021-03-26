//
//  constantforce.hpp
//  threege
//
//  Created by Maximilian Katzmann on 16.05.17.
//

#pragma once

#include <iostream>
#include <string>

#include "coordinate.h"
#include "force.h"

using namespace std;

class ConstantForce : public Force {
 public:
  /**
   * The constant force depends on the radius of the disk that the particles are
   * placed in.
   */
  double disk_radius;

  /**
   * Determines the maximum force that should be applied. I.e. the force that
   * is applied between the two vertices with the maximal distance.
   */
  double maximal_force;

  ConstantForce(double disk_radius, double maximal_force);

  ConstantForce();

  Coordinate force_from_coordinate_to_coordinate(const Coordinate& sender,
                                                 const Coordinate& receiver,
                                                 bool attractive);

  /**
   * The constant force uses the hyperbolic distance in the 3-dimensional
   * hyperbolic space.
   *
   * @param  coord1 A coordinate.
   * @param  coord2 Another coordinate.
   * @return        The distance between the coordinates in the 3-dimensional
   * hyperbolic space.
   */
  double distance_between_coordinates(const Coordinate& coord1,
                                      const Coordinate& coord2) const;

  virtual ~ConstantForce() = default;
  ConstantForce(const ConstantForce& other) = default;
  ConstantForce(ConstantForce&& other) = default;
  ConstantForce& operator=(const ConstantForce& other) = default;
  ConstantForce& operator=(ConstantForce&& other) = default;

  // private:
  //  double maximalDistanceBetweenParticles;
};

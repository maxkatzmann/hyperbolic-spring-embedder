//
//  force.hpp
//  threege
//
//  Created by Maximilian Katzmann on 09.05.17.
//
//  This class represent forces that act between vertices in the graph.
//  Note: This class is abstract. You need to derive from it in order
//  to define a force.
//

#pragma once

#include <iostream>
#include <string>

#include "coordinate.h"

using namespace std;

class Force {
 public:
  /**
   * The temperature may impact the force.
   */
  double temperature = 1.0;

  /**
   * Determines the factor that the temperature is multiplied by when decreasing
   * it after beeing applied. Should be a value in [0, 1].
   */
  double temperature_decrease = 0.9;

  Force();
  /**
   * The force acts from the sender on to the receiver. The coordinate that
   * is returned is the position of the receiver after the force has been
   * applied.
   * @param  sender     The particle propagating the force.
   * @param  receiver   The particle that the force acts on.
   * @param  attractive Determines whether the force is attractive (true) or
   * repulsive (false).
   * @return            The position of the receiver after the force has been
   * applied.
   */
  virtual Coordinate force_from_coordinate_to_coordinate(
      const Coordinate& sender, const Coordinate& receiver, bool attractive);

  virtual double radial_force_from_coordinate_to_coordinate(
      const Coordinate& sender, const Coordinate& receiver, bool attractive);
  /**
   * The distance between coordinates. May be used when applying forces.
   * @param  coord1 A coordinate.
   * @param  coord2 Another coordinate.
   * @return        The distance between the two coordinates.
   */
  virtual double distance_between_coordinates(const Coordinate& coord1,
                                              const Coordinate& coord2) const;
  /**
   * Determines the normal vector of a plane that is defined by the three passed
   * coordinates.
   * @param  coord1 A coordinate that lies on the plane.
   * @param  coord2 Another coordinate that lies on the plane.
   * @param  coord3 A third coordinate that lies on the plane.
   * @return        The unit vector that lies perpendicular on that plane.
   */

  static double hyperbolic_distance_between_coordinates(
      const Coordinate& coord1, const Coordinate& coord2);

  ~Force() = default;
  Force(const Force& other) = default;
  Force(Force&& other) = default;
  Force& operator=(const Force& other) = default;
  Force& operator=(Force&& other) = default;
};

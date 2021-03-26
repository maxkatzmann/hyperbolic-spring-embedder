//
//  coordinate.hpp
//  threege
//
//  Created by Maximilian Katzmann on 04.05.17.
//
//  This class represents a 3-dimensional cartesian coordinate
//  or a 3-dimensional vector.
//  It provides methods for adding coordinates, rotating them
//  around an axis (defined by another coordinate).

#pragma once

#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

class Coordinate {
 public:

  double x, y, z;

  Coordinate() : x(0.0), y(0.0), z(0.0) {};

  ~Coordinate() = default;
  Coordinate(const Coordinate& other) = default;
  Coordinate(Coordinate&& other) = default;

  Coordinate(double x, double y, double z) : x(x), y(y), z(z) {};

  double length() const;

  double phi() const;

  Coordinate cross_product(const Coordinate& other) const;

  double dot_product(const Coordinate& other) const;

  /**
   * Rotates the receiver around the passed
   * axis by the given angle according to the right-hand-rule.
   *
   * Rodrigues' formula is used to perform the rotation.
   *
   * Note: The vector defining the axis of rotation has to have unit length.
   *
   * @param  axis  A coordinate representing the UNIT vector that the orig
   * coordinate should be rotated about
   * @param  angle The angle in radians describing how far orig is rotated
   * around the axis.
   * @return       A coordinate representing the result of the rotation.
   */
  Coordinate coordinate_by_rotation(Coordinate& axis, double angle) const;

  /**
   * Converts the receiver which is expected to describe cartesian coordinates
   * to spherical coordinates.
   * @return A coordinate describing the receiver using spherical coordinates.
   */
  Coordinate cartesian_to_spherical() const;

  /**
   * Converts the receiver which expected to describe a spherical coordinate to
   * a cartesian coordinate.
   * @return A coordinate describing the receiver using cartesian coordinates.
   */
  Coordinate spherical_to_cartesian() const;

  /**
   * Determines the central angle between the receiver and other.
   * @note This method assumes that the receiver and other use spherical coordinates.
   *
   * @param other The spherical coordinate describing the point to which the central angle is to be obtained.
   * @return The central angle between the receiver and other in radians.
   */
  double central_angle_to_coordinate(const Coordinate& other) const;

  static Coordinate normal_vector_for_plane_from_points(const Coordinate &coord1,
                                                        const Coordinate &coord2,
                                                        const Coordinate &coord3);

  /**
   * Determines a random point in the 3-dimensional euclidean unit ball.
   *
   * Note: The point is not chosen uniformly. The probability of the point
   * lying on the surface of the sphere is ~48%.
   *
   * @param  seed Used for random generation.
   * @return A coordinate resembling the random point.
   */
  static Coordinate random_point_in_unit_ball(int seed);

  /**
   * Determines are random point on the surface of a 3-dimensional euclidean
   * unit ball.
   *
   * @param  seed Used for random generation.
   * @return      A 3-dimensional coordinate representing a random point on the
   * surface of the euclidean unit ball.
   */
  static Coordinate random_point_on_surface_of_unit_ball(int seed);

  std::string description() const {
    std::ostringstream oss;
    oss << std::fixed << std::showpoint;
    oss << std::setprecision(7);
    oss << "(" << x << ", " << y << ", " << z << ")";
    return oss.str();
  }

  Coordinate& operator=(const Coordinate& otVher) = default;
  Coordinate& operator=(Coordinate&& other) = default;

  bool operator==(const Coordinate& other) const;

  Coordinate operator+(const Coordinate& other) const;
  Coordinate operator+=(const Coordinate& other);
  Coordinate operator-(const Coordinate& other) const;
  Coordinate operator*(const double factor) const;
  Coordinate operator*=(const double factor);
  Coordinate operator/(const double facotr) const;
  Coordinate operator/=(const double factor);

};

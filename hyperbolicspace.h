//
//  hyperbolicspace.hpp
//  threege
//
//  Created by Maximilian Katzmann on 24.10.16.
//  Copyright © 2016 Max Katzmann. All rights reserved.
//

#ifndef hyperbolicspace_hpp
#define hyperbolicspace_hpp

#include <math.h>
#include <stdio.h>

#include <cmath>

#include "coordinate.h"
#include "force.h"

using namespace std;

class HyperbolicSpace {
 public:
  static double distance_between_native_coordinate(const Coordinate &coord1,
                                                   const Coordinate &coord2) {
    /**
     * Check whether one of the coordinates is the origin. If it is, we can
     * simply return the radius of the second coordinate.
     */
    Coordinate origin = Coordinate();
    if (coord1 == origin) {
      return coord2.x;
    }
    if (coord2 == origin) {
      return coord1.x;
    }

    double lambda1 = coord1.y;
    double lambda2 = coord2.y;
    double phi1 = coord1.z;
    double phi2 = coord2.z;

    double deltaLambda = abs(lambda1 - lambda2);

    /**
     * The central angle is determined using the following function, which is
     * said to be rather robust against numerical difficulties.
     */
    double deltaSigma = atan2(
        sqrt(
            (cos(phi2) * sin(deltaLambda)) * (cos(phi2) * sin(deltaLambda)) +
            (cos(phi1) * sin(phi2) - sin(phi1) * cos(phi2) * cos(deltaLambda)) *
                (cos(phi1) * sin(phi2) -
                 sin(phi1) * cos(phi2) * cos(deltaLambda))),
        (sin(phi1) * sin(phi2) + cos(phi1) * cos(phi2) * cos(deltaLambda)));

    /**
     * Now we now the angular distance between the to particles. Next up: the
     * radial coordinates.
     */

    double r1 = coord1.x;
    double r2 = coord2.x;

    return acosh(cosh(r1) * cosh(r2) - sinh(r1) * sinh(r2) * cos(deltaSigma));
  }

  static double euclideanRadialDistanceToHyperbolicRadialDistance(double d) {
    return 2.0 * atanh(d);
  }

  static double hyperbolicRadialDistanceToEuclideanRadialDistance(double d) {
    return tanh(d / 2.0);
  }

  static Coordinate coordinate_mirrored_on_xz_plane(Coordinate coordinate) {
    return Coordinate(coordinate.x, -(coordinate.y - (2.0 * M_PI)),
                      coordinate.z);
  }
  static Coordinate coordinate_mirrored_on_xy_plane(Coordinate coordinate) {
    return Coordinate(coordinate.x, coordinate.y, -coordinate.z);
  }

  static Coordinate translate_point_along_x_axis_natively(
      Coordinate point, double hyperbolicDistance) {
    Coordinate original_point = point;

    Coordinate reference_point;
    if (hyperbolicDistance > 0.0) {
      reference_point = Coordinate(abs(hyperbolicDistance), M_PI, 0.0);
    } else if (hyperbolicDistance < 0.0) {
      reference_point = Coordinate(abs(hyperbolicDistance), 0.0, 0.0);
    } else {
      return point;
    }

    if (original_point.y > M_PI) {
      point = HyperbolicSpace::coordinate_mirrored_on_xz_plane(point);
    }

    if (original_point.z < 0) {
      point = HyperbolicSpace::coordinate_mirrored_on_xy_plane(point);
    }

    double d = abs(hyperbolicDistance);
    double r = point.x;
    double phi = point.z;

    double r_prime = HyperbolicSpace::distance_between_native_coordinate(
        point, reference_point);

    double z = asinh(sinh(r) * sin(phi));
    double x = acosh(cosh(r) / cosh(z));
    double x_prime = acosh(cosh(r_prime) / cosh(z));
    double phi_prime = asin(sinh(z) / sinh(r_prime));
    double lambda_prime =
        acos((cosh(x_prime) * cosh(d) - cosh(x)) / (sinh(x_prime) * sinh(d)));

    if (hyperbolicDistance < 0.0) {
      lambda_prime = M_PI - lambda_prime;
    }

    if (lambda_prime != lambda_prime) {
      lambda_prime = point.y;
    }

    if (phi_prime != phi_prime) {
      phi_prime = point.z;
    }

    Coordinate translated_coordinate =
        Coordinate(r_prime, lambda_prime, phi_prime);
    if (original_point.y > M_PI) {
      translated_coordinate = HyperbolicSpace::coordinate_mirrored_on_xz_plane(
          translated_coordinate);
    }

    if (original_point.z < 0) {
      translated_coordinate = HyperbolicSpace::coordinate_mirrored_on_xy_plane(
          translated_coordinate);
    }

    return translated_coordinate;
  }

  /**
   *  ASSUMES CARTESIAN COORDINATES!
   */
  static Coordinate translate_point_by_vector(const Coordinate &point,
                                              const Coordinate &vector) {
    if (fabs(point.x - vector.x) < 0.0000000001 &&
        fabs(point.y - vector.y) < 0.0000000001 &&
        fabs(point.z - vector.z) < 0.0000000001) {
      return Coordinate();
    }

    if (vector.y == 0.0 && vector.z == 0.0) {
      Coordinate spherical_point = point.cartesian_to_spherical();
      Coordinate translated_coordinate =
          HyperbolicSpace::translate_point_along_x_axis_natively(
              spherical_point, vector.x);
      return translated_coordinate.spherical_to_cartesian();
    }

    /**
     *  First we rotate everything such that the vector lies on the x-axis.
     */
    Coordinate reference_point = Coordinate(vector.length(), 0.0, 0.0);

    Coordinate origin = Coordinate();
    Coordinate rotation_axis = Coordinate::normal_vector_for_plane_from_points(
        origin, vector, reference_point);

    Coordinate spherical_vector = vector.cartesian_to_spherical();
    Coordinate spherical_reference_point =
        reference_point.cartesian_to_spherical();

    double rotation_angle =
        spherical_vector.central_angle_to_coordinate(spherical_reference_point);

    Coordinate rotated_point =
        point.coordinate_by_rotation(rotation_axis, rotation_angle);

    /**
     *  Now we translate everything such that the vector is placed on the
     * origin.
     */
    Coordinate rotated_point_spherical = rotated_point.cartesian_to_spherical();
    Coordinate translated_point_spherical =
        HyperbolicSpace::translate_point_along_x_axis_natively(
            rotated_point_spherical, -reference_point.x);

    Coordinate translated_point =
        translated_point_spherical.spherical_to_cartesian();

    if (translated_point.x != translated_point.x) {
      std::cout << "Here!" << std::endl;
    }

    /**
     *  Finally we rotate everything back.
     */
    Coordinate inverted_rotation_axis = rotation_axis * -1;
    return translated_point.coordinate_by_rotation(inverted_rotation_axis,
                                                   rotation_angle);
  }

  /**
   * Determines the weighted centroid on the line between the two passed
   * coordinates.
   *
   * @note It is assummed that the coordinates are present in spherical
   * coordinates.
   *
   * @param coord1 The spherical representation of the first coordinate.
   * @param weight1 The weight of the first coordinate.
   * @param coord2 The spherical representation of the second coordinate.
   * @param weight2 The weight of the second coordinate.
   * @return The weighted centroid on the line between the two coordinates.
   */
  static Coordinate weighted_centroid_between_coordinates(
      const Coordinate &coord1, double weight1, const Coordinate &coord2,
      double weight2) {
    if (coord1.x == 0.0) {
      double new_radius = coord2.x * (weight2 / (weight1 + weight2));
      return Coordinate(new_radius, coord2.y, coord2.z);
    }

    if (coord2.x == 0.0) {
      double new_radius = coord1.x * (weight1 / (weight1 + weight2));
      return Coordinate(new_radius, coord1.y, coord1.z);
    }

    double central_angle = coord1.central_angle_to_coordinate(coord2);

    double difference_to_M_PI = fabs(central_angle - M_PI);
    if (difference_to_M_PI < 0.0000000001) {
      double distance = coord1.x + coord2.x;
      double weighted_distance_from_q =
          distance * (weight1 / (weight1 + weight2));

      double new_radius = coord2.x - weighted_distance_from_q;

      if (new_radius < 0.0) {
        return Coordinate(fabs(new_radius), coord1.y, coord1.z);
      } else {
        return Coordinate(new_radius, coord2.y, coord2.z);
      }
    }

    if (fabs(central_angle) < 0.0000000001) {
      double distance = fabs(coord1.x - coord2.x);
      double weighted_distance_from_q =
          distance * (weight1 / (weight1 + weight2));

      if (coord2.x > coord1.x) {
        return Coordinate(coord2.x - weighted_distance_from_q, coord2.y,
                          coord2.z);
      } else {
        return Coordinate(coord2.x + weighted_distance_from_q, coord2.y,
                          coord2.z);
      }
    }

    double d =
        HyperbolicSpace::distance_between_native_coordinate(coord1, coord2);
    double d_q = d * (weight1 / (weight1 + weight2));
    double r_p = coord1.x;
    double r_q = coord2.x;

    double gamma =
        acos((cosh(d) * cosh(r_q) - cosh(r_p)) / (sinh(d) * sinh(r_q)));

    double r_c =
        acosh(cosh(d_q) * cosh(r_q) - sinh(d_q) * sinh(r_q) * cos(gamma));
    double alpha =
        acos((cosh(r_c) * cosh(r_q) - cosh(d_q)) / (sinh(r_c) * sinh(r_q)));

    /**
     *  Now we rotate q towards p by alpha and afterwards set the radius of
     *  the resulting point to r_c.
     */
    Coordinate coord_p_cartesian = coord1.spherical_to_cartesian();
    Coordinate coord_q_cartesian = coord2.spherical_to_cartesian();
    Coordinate origin = Coordinate();
    Coordinate rotation_axis = Coordinate::normal_vector_for_plane_from_points(
        origin, coord_q_cartesian, coord_p_cartesian);

    Coordinate rotated_q =
        coord_q_cartesian.coordinate_by_rotation(rotation_axis, alpha);
    Coordinate rotated_q_spherical = rotated_q.cartesian_to_spherical();
    rotated_q_spherical.x = r_c;

    return rotated_q_spherical;
  }
};

#endif /* hyperbolicspace_hpp */

//
//  embedding.hpp
//  threege
//
//  Created by Maximilian Katzmann on 04.05.17.
//

#pragma once

#include <iostream>
#include <string>

#include "coordinate.h"
#include "force.h"
#include "graph.h"

class Embedding {
 public:
  /**
   * The graph that this embedding represents.
   */
  Graph G;

  /**
   * The vector containing the coordinate for each vertex.
   */
  std::vector<Coordinate> points;

  /**
   * The velocity of the vertices which determines how much of the
   * previous forces are again part of the current force.
   */
  double velocity = 0.5;

  /**
   * Vectors storing information about forces during the embedding process.
   */
  std::vector<double> potential_history;
  std::vector<double> temperate_force_history;
  std::vector<double> average_plane_distance_history;
  std::vector<Coordinate> current_forces;
  std::vector<double> radial_potential_history;
  std::vector<double> temperate_radial_potential_history;

  /**
   * Constructors
   */
  Embedding() : G(0){};

  Embedding(const Graph& G, std::vector<Coordinate>& points);

  ~Embedding() = default;
  Embedding(const Embedding& other) = default;

  void read_embedding_from_file(
      const std::string& filename,
      std::unordered_map<std::string, int>& label_to_node);

  void read_embedding_from_file_3D(
      const std::string& filename,
      std::unordered_map<std::string, int>& label_to_node);

  double max_radius() const;

  void apply_permutation(const std::vector<int>& permutation);

  /**
   * Determines the force that acts on the passed vertex.
   * @param  force  The force that is being applied.
   * @param  vertex The vertex that the force acts on.
   * @return        A vector representing the force acts on the vertex.
   */
  Coordinate combined_force_on_vertex(Force& force, int vertex);

  /**
   * Determines the force that acts on the passed vertex radially.
   * @param  force  The force that is being applied.
   * @param  vertex The vertex that the force acts on.
   * @return        A vector representing the force acts on the vertex.
   */
  double combined_radial_force_on_vertex(Force& force, int vertex);

  /**
   * Determines an approximation of the force that acts on the passed vertex.
   * @param  force                 The force that is being applied.
   * @param  vertex                The vertex that the force acts on.
   * @param pooling_threshold If a cells distance to the vertex is larger
   * than this value, the cells center will be used instead of the coordinates
   * of the vertices in the cell.
   * @param distance_threshhold_2 If a cells distance to the vertex is larger
   * than this value, the children of the cell will not be processed at all.
   * @return                       A vector representing the force acts on the
   * vertex.
   */
  Coordinate combined_approximate_force_on_vertex(Force& force, int vertex);

  /**
   * Determines an approximation of the force that acts on the passed vertex
   * radially.
   * @param  force  The force that is being applied.
   * @param  vertex The vertex that the force acts on.
   * @param pooling_threshold If a cells distance to the vertex is larger
   * than this value, the cells center will be used instead of the coordinates
   * of the vertices in the cell.
   * @param distance_threshhold_2 If a cells distance to the vertex is larger
   * than this value, the children of the cell will not be processed at all.
   * @return        A double denoting the change in radius.
   */
  double combined_approximate_radial_force_on_vertex(Force& force, int vertex);

  /**
   * Applies the passed force between all pairs of vertices.
   * @param force A force object that determines which forces are applied.
   */
  void apply_force_to_vertices(Force& force);

  void apply_sampled_force_to_vertices(
      Force& force, const double sample_balancing_coefficient);

  void apply_force_to_vertices_and_flatten_towards_plane_with_normal(
      Force& force, Force& plane_force, Coordinate& unit_plane_normal);

  /**
   * Applies the passed force between vertices using a geometric data structure
   * for speed up.
   * @param force                 A force object that determines which forces
   * are applied.
   * @param plane_Force            The force the pulls the vertices towards the
   * plane.
   * @param planeNormal           The normal vector representing the plane that
   * the vertices should be pulled towards.
   * @param pooling_threshold If a cells distance to the vertex is larger
   * than this value, the cells center will be used instead of the coordinates
   * of the vertices in the cell.
   * @param distance_threshhold_2 If a cells distance to the vertex is larger
   * than this value, the children of the cell will not be processed at all.
   */
  void
  apply_approximate_force_to_vertices_and_flatten_towards_plane_with_normal(
      Force& force, Force& plane_force, Coordinate& unit_plane_normal);

  void apply_radial_force_to_vertices(Force& force);

  /**
   * Applies forces that only act on the radii of the vertices.
   * @param force                 The force to act on the radii.
   * @param pooling_threshold If a cells distance to the vertex is larger
   * than this value, the cells center will be used instead of the coordinates
   * of the vertices in the cell.
   * @param distance_threshhold_2 If a cells distance to the vertex is larger
   * than this value, the children of the cell will not be processed at all.
   */
  void apply_approximate_radial_force_to_vertices(Force& force);

  /**
   * Determines the expected degree of the current node depending on
   * the radii of all other nodes under the assumption that all nodes
   * have an angle that is uniformly chosen.
   */
  double current_expected_degree_of_vertex(int vertex);

  Coordinate fitting_plane() const;

  Coordinate euclidean_centroid() const;

  Coordinate approximate_centroid() const;

  void center();

  void move_vector_to_origin(const Coordinate& vector);

  double distance_from_point_to_plane_with_normal(Coordinate& point,
                                                  Coordinate& normal);

  double average_distance_of_embedding_to_plane_with_normal(Coordinate& normal);

  void convert_to_2D_embedding();

  bool is_stable() const;
  bool is_stable_with_threshold(double threshold) const;
  bool is_radially_stable() const;

  void draw(const char* filename, const vector<int>& nodes,
            int nodeDegreeThreshhold, const string& add_to_svg, int height,
            int width) const;

  std::string edge_histogram();

  void edge_length_histogram(std::vector<std::pair<int, int>>& histogram);

  void get_phi_phi_plot(vector<double>& other_phi_values,
                        vector<pair<double, double>>& phi_phi_plot);
  void get_r_r_plot(vector<double>& other_r_values,
                    vector<pair<double, double>>& r_r_plot);

  void get_phi_values(vector<double>& phi_values);
  void get_r_values(vector<double>& r_values);

  double log_likelihood() const;

  string to_string(bool useLabels = true) const;
  string to_string_spherical(bool useLabels = true) const;
  string to_string_2D(bool useLabels = true) const;

  Embedding(Embedding&& other) = default;
  Embedding& operator=(const Embedding& other) = default;
  Embedding& operator=(Embedding&& other) = default;
};

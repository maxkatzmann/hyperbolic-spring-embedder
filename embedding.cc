#include "embedding.h"

#include <glog/logging.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

#include <Eigen/Eigen>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <queue>
#include <random>
#include <sstream>

#include "hyperbolicspace.h"

using namespace Eigen;

#define FOR(i, n) for (int(i) = 0; (i) < (n); (i)++)

Embedding::Embedding(const Graph& G, vector<Coordinate>& points)
    : G(G), points(points) {
  current_forces = vector<Coordinate>(G.n, Coordinate());
}

double Embedding::max_radius() const {
  double max_radius = 0.0;

  FOR(i, G.n) {
    double radius = points[i].length();
    max_radius = std::max(radius, max_radius);
  }

  return max_radius;
}

void Embedding::apply_permutation(const std::vector<int>& permutation) {
  std::vector<Coordinate> new_points(G.n);

  FOR(i, G.n) { new_points[i] = points[permutation[i]]; }
  points = new_points;
}

std::string Embedding::to_string(bool use_labels) const {
  std::ostringstream oss;
  FOR(i, G.n) {
    if (use_labels) {
      oss << G.labels[i] << "\t" << points[i].x << "\t" << points[i].y << "\t"
          << points[i].z << std::endl;
    } else {
      oss << i << "\t" << points[i].x << "\t" << points[i].y << "\t"
          << points[i].z << std::endl;
    }
  }

  return oss.str();
}

string Embedding::to_string_spherical(bool use_labels) const {
  ostringstream oss;
  FOR(i, G.n) {
    Coordinate spherical_point = points[i].cartesian_to_spherical();
    if (use_labels) {
      oss << G.labels[i] << "\t(" << spherical_point.y * 180.0 / M_PI << ", "
          << spherical_point.z * 180.0 / M_PI << ")" << std::endl;
    } else {
      oss << i << "\t(" << spherical_point.y * 180.0 / M_PI << ", "
          << spherical_point.z * 180.0 / M_PI << ")" << std::endl;
    }
  }

  return oss.str();
}

string Embedding::to_string_2D(bool use_labels) const {
  ostringstream oss;
  FOR(i, G.n) {
    const auto point = points[i];
    const double radius = point.length();
    const double angle = point.phi();
    const double angleDegrees = angle / M_PI * 180.0;

    if (use_labels) {
      oss << G.labels[i] << "\t" << radius << "\t" << angleDegrees << "\n";
    } else {
      oss << i << "\t" << radius << "\t" << angleDegrees << "\n";
    }
  }

  return oss.str();
}

Coordinate Embedding::combined_force_on_vertex(Force& force, int vertex) {
  /**
   * Collect all the forces that act on the vertex.
   */
  Coordinate combinedForce = current_forces[vertex] * 0.25;
  FOR(j, G.n) {
    if (vertex != j) {
      /**
       * We now determine the force that vertex i applied on j. That is,
       * i is the sender and j is the receiver.
       */
      bool attractive = G.adjacent(vertex, j);
      Coordinate iCoord = points[vertex];
      Coordinate jCoord = points[j];

      Coordinate force_vector =
          force.force_from_coordinate_to_coordinate(jCoord, iCoord, attractive);
      combinedForce += force_vector;
    }
  }

  return combinedForce;
}

Coordinate Embedding::combined_approximate_force_on_vertex(Force& force,
                                                           int vertex) {
  // To simulate velocity, we consider the current_force again.
  Coordinate combined_force = current_forces[vertex] * this->velocity;

  return combined_force;
}

double Embedding::combined_radial_force_on_vertex(Force& force, int vertex) {
  /**
   * The force is determined as follows. If the currently expected
   * degree of the vertex (only depending on the radii of the other
   * vertices) is too small, we move the vertex inwards. If this expected
   * degree is too large, we move it outwards.
   */
  double currently_expected_degree = current_expected_degree_of_vertex(vertex);
  int actual_degree = G.edges[vertex].size();

  if (currently_expected_degree < actual_degree) {
    return -force.temperature;
  } else if (currently_expected_degree > actual_degree) {
    return force.temperature;
  } else {
    return 0.0;
  }
}

double Embedding::combined_approximate_radial_force_on_vertex(Force& force,
                                                              int vertex) {
  double combined_force = 0.0;

  return combined_force;
}

void Embedding::apply_force_to_vertices(Force& force) {
  int n = G.n;
  /**
   * We first determine the new positioins for all vertices and save them
   * afterwards. Else the new positions will already affect the current
   * iteration.
   */
  vector<Coordinate> new_coordinate_for_vertex(G.n);

  double temperate_force_sum = 0.0;
  FOR(i, n) {
    Coordinate combined_force = combined_force_on_vertex(force, i);
    current_forces[i] = combined_force;
    Coordinate rotation_axis = combined_force / combined_force.length();
    double rotation_angle = combined_force.length();

    temperate_force_sum += rotation_angle;

    new_coordinate_for_vertex[i] =
        points[i].coordinate_by_rotation(rotation_axis, rotation_angle);
  }

  FOR(i, n) { points[i] = new_coordinate_for_vertex[i]; }

  potential_history.push_back(temperate_force_sum / force.temperature);
  temperate_force_history.push_back(temperate_force_sum);
}

void Embedding::apply_sampled_force_to_vertices(
    Force& force, const double sample_balancing_coefficient) {
  double temperate_force_sum = 0.0;

  /**
   * We first determine the new positioins for all vertices and save them
   * afterwards. Else the new positions will already affect the current
   * iteration.
   */
  vector<Coordinate> new_coordinate_for_vertex(G.n);

  FOR(i, G.n) {
    Coordinate combined_force = combined_approximate_force_on_vertex(force, i);

    current_forces[i] = combined_force;

    double rotation_angle = combined_force.length();
    Coordinate rotation_axis = combined_force / rotation_angle;

    temperate_force_sum += rotation_angle;

    /**
     * Now we obtain the new coordinate by applying the rotation.
     */
    new_coordinate_for_vertex[i] =
        points[i].coordinate_by_rotation(rotation_axis, rotation_angle);
  }

  /**
   * Updating the cell data structure
   */
  FOR(i, G.n) { points[i] = new_coordinate_for_vertex[i]; }

  potential_history.push_back(temperate_force_sum / force.temperature);
  temperate_force_history.push_back(temperate_force_sum);
}

void Embedding::apply_force_to_vertices_and_flatten_towards_plane_with_normal(
    Force& force, Force& plane_force, Coordinate& unit_plane_normal) {
  int n = G.n;

  /**
   * We first determine the new positioins for all vertices and save them
   * afterwards. Else the new positions will already affect the current
   * iteration.
   */
  vector<Coordinate> new_coordinate_for_vertex(G.n);

  double temperate_force_sum = 0.0;

  double max_radius = this->max_radius();

  FOR(i, n) {
    Coordinate combined_force = combined_force_on_vertex(force, i);

    current_forces[i] = combined_force;
    double radius = points[i].length();

    /**
     *  Now we add the force that pulls the vertex towards the plane.
     *  Therefore, we get the coordinate of the vertex on the plane.
     *
     *  Meaning we take the radius of the vertex and project its angular
     *  coordinate onto the plane.
     */
    Coordinate rotation_axis_onto_plane =
        Coordinate::normal_vector_for_plane_from_points(
            Coordinate(), unit_plane_normal, points[i]);

    rotation_axis_onto_plane /= rotation_axis_onto_plane.length();
    Coordinate point_on_normal_vector = unit_plane_normal * radius;

    Coordinate point_on_plane = point_on_normal_vector.coordinate_by_rotation(
        rotation_axis_onto_plane, M_PI_2);
    double temp_temperature = plane_force.temperature;
    plane_force.temperature *= (radius / max_radius);
    Coordinate plane_attraction_force =
        plane_force.force_from_coordinate_to_coordinate(point_on_plane,
                                                        points[i], true);
    plane_force.temperature = temp_temperature;

    combined_force += plane_attraction_force;

    Coordinate rotation_axis = combined_force / combined_force.length();
    double rotation_angle = combined_force.length();
    temperate_force_sum += plane_attraction_force.length();

    new_coordinate_for_vertex[i] =
        points[i].coordinate_by_rotation(rotation_axis, rotation_angle);
  }

  FOR(i, n) { points[i] = new_coordinate_for_vertex[i]; }

  potential_history.push_back(temperate_force_sum / force.temperature);
  temperate_force_history.push_back(temperate_force_sum);
  average_plane_distance_history.push_back(
      average_distance_of_embedding_to_plane_with_normal(unit_plane_normal));
}

void Embedding::
    apply_approximate_force_to_vertices_and_flatten_towards_plane_with_normal(
        Force& force, Force& plane_force, Coordinate& unit_plane_normal) {
  /**
   * We first determine the new positioins for all vertices and save them
   * afterwards. Else the new positions will already affect the current
   * iteration.
   */
  vector<Coordinate> new_coordinate_for_vertex(G.n);

  double temperate_force_sum = 0.0;

  double max_radius = this->max_radius();

  FOR(i, G.n) {
    Coordinate combined_force = combined_approximate_force_on_vertex(force, i);

    current_forces[i] = combined_force;
    double radius = points[i].length();

    /**
     *  Now we add the force that pulls the vertex towards the plane.
     *  Therefore, we get the coordinate of the vertex on the plane.
     *
     *  Meaning we take the radius of the vertex and project its angular
     *  coordinate onto the plane.
     */
    Coordinate rotation_axis_onto_plane =
        Coordinate::normal_vector_for_plane_from_points(
            Coordinate(), unit_plane_normal, points[i]);

    rotation_axis_onto_plane /= rotation_axis_onto_plane.length();
    Coordinate point_on_normal_vector = unit_plane_normal * radius;

    Coordinate point_on_plane = point_on_normal_vector.coordinate_by_rotation(
        rotation_axis_onto_plane, M_PI_2);
    double temp_temperature = plane_force.temperature;
    plane_force.temperature *= (radius / max_radius);
    Coordinate plane_attraction_force =
        plane_force.force_from_coordinate_to_coordinate(point_on_plane,
                                                        points[i], true);
    plane_force.temperature = temp_temperature;

    combined_force += plane_attraction_force;

    Coordinate rotation_axis = combined_force / combined_force.length();
    double rotation_angle = combined_force.length();

    temperate_force_sum += plane_attraction_force.length();

    new_coordinate_for_vertex[i] =
        points[i].coordinate_by_rotation(rotation_axis, rotation_angle);
  }

  FOR(i, G.n) { points[i] = new_coordinate_for_vertex[i]; }

  potential_history.push_back(temperate_force_sum / force.temperature);
  temperate_force_history.push_back(temperate_force_sum);
  average_plane_distance_history.push_back(
      average_distance_of_embedding_to_plane_with_normal(unit_plane_normal));
}

void Embedding::apply_radial_force_to_vertices(Force& force) {
  int n = G.n;

  /**
   * We first determine the new positioins for all vertices and save them
   * afterwards. Else the new positions will already affect the current
   * iteration.
   */
  vector<Coordinate> new_coordinate_for_vertex(G.n);

  double radial_potential_sum = 0.0;
  FOR(i, n) {
    double combined_force = combined_radial_force_on_vertex(force, i);

    radial_potential_sum += fabs(combined_force);

    double current_radius = points[i].length();
    double new_radius = current_radius + combined_force;
    if (new_radius < 0.0) {
      new_radius = 0.0;
    }

    Coordinate current_point = points[i] / current_radius;
    new_coordinate_for_vertex[i] = current_point * new_radius;
  }

  FOR(i, G.n) {
    /**
     * Updating the cell data structure.
     */
    // Coordinate current_spherical = points[i].cartesian_to_spherical();
    // Coordinate new_spherical =
    //     new_coordinate_for_vertex[i].cartesian_to_spherical();
    // cell_structure.remove_node_with_position(i, current_spherical);
    // int new_cell = cell_structure.add_node_with_position(i, new_spherical);
    // cell_for_vertex[i] = new_cell;
    points[i] = new_coordinate_for_vertex[i];
  }

  radial_potential_history.push_back(radial_potential_sum / force.temperature);
  temperate_radial_potential_history.push_back(radial_potential_sum);
}

void Embedding::apply_approximate_radial_force_to_vertices(Force& force) {
  /**
   * We first determine the new positioins for all vertices and save them
   * afterwards. Else the new positions will already affect the current
   * iteration.
   */
  vector<Coordinate> new_coordinate_for_vertex(G.n);

  double radial_potential_sum = 0.0;
  FOR(i, G.n) {
    double combined_force =
        combined_approximate_radial_force_on_vertex(force, i);

    radial_potential_sum += fabs(combined_force);

    double current_radius = points[i].length();
    double new_radius = current_radius + combined_force;
    if (new_radius > 0.0 && current_radius > 0.0) {
      new_coordinate_for_vertex[i] = (points[i] / current_radius) * new_radius;
    }
  }

  FOR(i, G.n) {
    /**
     * Updating the cell data structure.
     */
    points[i] = new_coordinate_for_vertex[i];
  }

  radial_potential_history.push_back(radial_potential_sum / force.temperature);
  temperate_radial_potential_history.push_back(radial_potential_sum);
}

void Embedding::read_embedding_from_file(
    const string& filename, unordered_map<std::string, int>& label_to_node) {
  std::ifstream file;
  file.open(filename.c_str());
  CHECK(file.good()) << "Maybe the file " << filename << " doesn't exist?";

  points.resize(G.n);

  // read all nodes
  string current_label, current_x, current_y;
  while (file >> current_label && file >> current_x && file >> current_y) {
    if ((file >> std::ws).peek() == '#') {
      // ignore rest of the line
      std::getline(file, current_label);
      continue;
    }

    current_label = current_label.substr(0, current_label.find("/"));
    current_x = current_x.substr(0, current_x.find("/"));
    // current_x.erase(std::remove(current_x.begin(), current_x.end(), '('),
    //                 current_x.end());
    current_y = current_y.substr(0, current_y.find("/"));
    if ((label_to_node).find(current_label) != label_to_node.end()) {
      double x = atof(current_x.c_str());
      double y = atof(current_y.c_str());
      points[label_to_node[current_label]] = Coordinate(x, y, 0.0);
    }
  }
  file.close();
}

void Embedding::read_embedding_from_file_3D(
    const string& filename, unordered_map<std::string, int>& label_to_node) {
  std::ifstream file;
  file.open(filename.c_str());
  CHECK(file.good()) << "Maybe the file " << filename << " doesn't exist?";

  points.resize(G.n);

  // read all nodes
  string current_label, current_x, current_y, current_z;
  while (file >> current_label && file >> current_x && file >> current_y &&
         file >> current_z) {
    if ((file >> std::ws).peek() == '#') {
      // ignore rest of the line
      std::getline(file, current_label);
      continue;
    }

    current_label = current_label.substr(0, current_label.find("/"));
    current_x = current_x.substr(0, current_x.find("/"));
    // current_x.erase(std::remove(current_x.begin(), current_x.end(), '('),
    //                 current_x.end());
    current_y = current_y.substr(0, current_y.find("/"));
    current_z = current_z.substr(0, current_z.find("/"));
    if ((label_to_node).find(current_label) != label_to_node.end()) {
      double x = atof(current_x.c_str());
      double y = atof(current_y.c_str());
      double z = atof(current_z.c_str());
      points[label_to_node[current_label]] = Coordinate(x, y, z);
    }
  }
  file.close();
}

void Embedding::draw(const char* filename, const vector<int>& nodes,
                     int nodeDegreeThreshhold, const string& add_to_svg,
                     int height, int width) const {
  int add = 0;

  FILE* file = fopen(filename, "w");
  if (file == NULL) {
    LOG(WARNING) << "Could not open file. " << strerror(errno);
    LOG(WARNING) << "Current working dir=" << getcwd(NULL, 0);
    ;
  }

  vector<std::string> colors = {"red",
                                "aliceblue",
                                "antiquewhite",
                                "aqua",
                                "aquamarine",
                                "azure",
                                "beige",
                                "bisque",
                                "blanchedalmond",
                                "blue",
                                "blueviolet",
                                "brown",
                                "burlywood",
                                "cadetblue",
                                "chartreuse",
                                "chocolate",
                                "coral",
                                "cornflowerblue",
                                "cornsilk",
                                "crimson",
                                "cyan",
                                "darkblue",
                                "darkcyan",
                                "darkgoldenrod",
                                "darkgray",
                                "darkgreen",
                                "darkgrey",
                                "darkkhaki",
                                "darkmagenta",
                                "darkolivegreen",
                                "darkorange",
                                "darkorchid",
                                "darkred",
                                "darksalmon",
                                "darkseagreen",
                                "darkslateblue",
                                "darkslategray",
                                "darkslategrey",
                                "darkturquoise",
                                "darkviolet",
                                "deeppink",
                                "deepskyblue",
                                "dimgray",
                                "dimgrey",
                                "dodgerblue",
                                "firebrick",
                                "floralwhite",
                                "forestgreen",
                                "fuchsia",
                                "gainsboro",
                                "ghostwhite",
                                "gold",
                                "goldenrod",
                                "gray",
                                "green",
                                "greenyellow",
                                "grey",
                                "honeydew",
                                "hotpink",
                                "indianred",
                                "indigo",
                                "ivory",
                                "khaki",
                                "lavender",
                                "lavenderblush",
                                "lawngreen",
                                "lemonchiffon",
                                "lightblue",
                                "lightcoral",
                                "lightcyan",
                                "lightgoldenrodyellow",
                                "lightgray",
                                "lightgreen",
                                "lightgrey",
                                "lightpink",
                                "lightsalmon",
                                "lightseagreen",
                                "lightskyblue",
                                "lightslategray",
                                "lightslategrey",
                                "lightsteelblue",
                                "lightyellow",
                                "lime",
                                "limegreen",
                                "linen",
                                "magenta(Safe",
                                "16=fuchsia",
                                "Hex3)",
                                "maroon",
                                "mediumaquamarine",
                                "mediumblue",
                                "mediumorchid",
                                "mediumpurple",
                                "mediumseagreen",
                                "mediumslateblue",
                                "mediumspringgreen",
                                "mediumturquoise",
                                "mediumvioletred",
                                "midnightblue",
                                "mintcream",
                                "mistyrose",
                                "moccasin",
                                "navajowhite",
                                "navy",
                                "oldlace",
                                "olive",
                                "olivedrab",
                                "orange",
                                "orangered",
                                "orchid",
                                "palegoldenrod",
                                "palegreen",
                                "paleturquoise",
                                "palevioletred",
                                "papayawhip",
                                "peachpuff",
                                "peru",
                                "pink",
                                "plum",
                                "powderblue",
                                "purple",
                                "rosybrown",
                                "royalblue",
                                "saddlebrown",
                                "salmon",
                                "sandybrown",
                                "seagreen",
                                "seashell",
                                "sienna",
                                "silver",
                                "skyblue",
                                "slateblue",
                                "slategray",
                                "slategrey",
                                "snow",
                                "springgreen",
                                "steelblue",
                                "tan",
                                "teal",
                                "thistle",
                                "tomato",
                                "turquoise",
                                "violet",
                                "wheat",
                                "white",
                                "whitesmoke",
                                "yellow",
                                "yellowgreen"};

  fprintf(file,
          "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n<!DOCTYPE svg PUBLIC "
          "\"-//W3C//DTD SVG 1.1//EN\" "
          "\"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd\">\n\n<svg "
          "xmlns=\"http://www.w3.org/2000/svg\"\nxmlns:xlink=\"http://"
          "www.w3.org/1999/xlink\" "
          "xmlns:ev=\"http://www.w3.org/2001/xml-events\"\nversion=\"1.1\" "
          "baseProfile=\"full\"\nwidth=\"%d\" height=\"%d\">\n\n",
          width + 2 * add, height + 2 * add);

  // comment in/out to see edges
  int num_nodes = 0;
  double R = this->max_radius();
  FOR(i, G.n)
  if (nodes[i] > 0) ++num_nodes;
  if (num_nodes <= 500000) {
    FOR(i, G.n) {
      FOR(jj, (int)G.edges[i].size()) {
        int j = G.edges[i][jj];
        if (!nodes[i] || !nodes[j]) continue;

        // skip long edges
        // TODO remove
        // if (HYPT::dist(pts[i], pts[j]) > 1.5 * R) continue;
        //
        // double opacity =
        //     std::min(1.0, (1.5 - Force::hyperbolicDistanceBetweenCoordinates(
        //                              points[i], points[j]) /
        //                    (R * 2.0)));
        double opacity = 1.0;

        bool doubleedge = false;
        FOR(k, jj)
        if (G.edges[i][k] == j) doubleedge = true;
        double x1 =
            add +
            0.5 * width * (1. + cos(points[i].phi()) * points[i].length() / R);
        double x2 =
            add +
            0.5 * width * (1. + cos(points[j].phi()) * points[j].length() / R);
        double y1 =
            add +
            0.5 * height * (1. + sin(points[i].phi()) * points[i].length() / R);
        double y2 =
            add +
            0.5 * height * (1. + sin(points[j].phi()) * points[j].length() / R);
        fprintf(file,
                "<line x1=\"%.8f\" y1=\"%.8f\" x2=\"%.8f\" y2=\"%.8f\" "
                "stroke=\"%s\" stroke-width=\"0.3\" opacity=\"%.8f\"/>\n",
                x1, y1, x2, y2, doubleedge ? "blue" : "gray", opacity);
      }
    }
  }

  FOR(i, G.n) {
    CHECK(nodes[i] < colors.size() + 2);
    if (nodes[i] == 0) continue;
    if (nodes[i] == 1)
      fprintf(file,
              "<circle cx=\"%.8f\" cy=\"%.8f\" r=\"%.8f\" fill=\"black\" "
              "stroke=\"black\" stroke-width=\"0\"/>\n",
              add + 0.5 * width *
                        (1. + cos(points[i].phi()) * points[i].length() / R),
              add + 0.5 * height *
                        (1. + sin(points[i].phi()) * points[i].length() / R),
              0.0015 * min(width, height));
    else
      fprintf(file,
              "<circle cx=\"%.8f\" cy=\"%.8f\" r=\"%.8f\" fill=\"%s\" "
              "stroke=\"%s\" stroke-width=\"0\"/>\n",
              add + 0.5 * width *
                        (1. + cos(points[i].phi()) * points[i].length() / R),
              add + 0.5 * height *
                        (1. + sin(points[i].phi()) * points[i].length() / R),
              0.003 * min(width, height), colors[nodes[i] - 2].c_str(),
              colors[nodes[i] - 2].c_str());
    //    if (num_nodes < 5000)

    if (G.edges[i].size() > nodeDegreeThreshhold) {
      fprintf(file,
              "<text x=\"%.8f\" y=\"%.8f\" font-family=\"sans-serif\" "
              "font-size=\"%dpx\" fill=\"red\" >%s</text>\n",
              add +
                  0.5 * width *
                      (1. + cos(points[i].phi()) * points[i].length() / R) +
                  0.002 * min(width, height),
              add + 0.5 * height *
                        (1. + sin(points[i].phi()) * points[i].length() / R),
              15, G.labels[i].c_str());
    }
  }

  // Boundary of hypRG
  fprintf(file,
          "<circle cx=\"%.8f\" cy=\"%.8f\" r=\"%.8f\" style=\"stroke: blue; "
          "fill: none;\"/>\n",
          add + 0.5 * width, add + 0.5 * height, 0.5 * min(width, height));

  fwrite(add_to_svg.c_str(), sizeof(char), add_to_svg.size(), file);

  fprintf(file, "\n</svg>\n");
  fclose(file);
}

void Embedding::get_phi_phi_plot(vector<double>& other_phi_values,
                                 vector<pair<double, double>>& phi_phi_plot) {
  FOR(i, G.n) {
    double phi_value = points[i].phi();
    double other_phi_value = other_phi_values[i];

    phi_phi_plot[i] = pair<double, double>(other_phi_value, phi_value);
  }
}

void Embedding::get_r_r_plot(vector<double>& other_r_values,
                             vector<pair<double, double>>& r_r_plot) {
  FOR(i, G.n) {
    double r_value = points[i].length();
    double other_r_value = other_r_values[i];

    r_r_plot[i] = pair<double, double>(other_r_value, r_value);
  }
}

void Embedding::get_phi_values(vector<double>& phi_values) {
  FOR(i, G.n) { phi_values[i] = points[i].phi(); }
}

void Embedding::get_r_values(vector<double>& r_values) {
  FOR(i, G.n) { r_values[i] = points[i].length(); }
}

double Embedding::log_likelihood() const {
  double R = this->max_radius();
  double T = 0.1;
  double summed_log_likelihood = 0.0;
  FOR(i, G.n) {
    FOR(j, G.n) {
      if (i != j) {
        double distance = Force::hyperbolic_distance_between_coordinates(
            points[i], points[j]);
        if (G.adjacent(i, j)) {
          summed_log_likelihood +=
              log((1.0 / (1.0 + exp((1.0 / (2.0 * T)) * (distance - R)))));
        } else {
          if (distance < R) {
            double inverted_distance = 2.0 * R - distance;
            summed_log_likelihood += log(
                1.0 / (1.0 + exp((1.0 / (2.0 * T)) * (inverted_distance - R))));
          } else {
            summed_log_likelihood += log(
                1.0 - (1.0 / (1.0 + exp((1.0 / (2.0 * T)) * (distance - R)))));
          }
        }
      }
    }
  }

  return summed_log_likelihood;
}

string Embedding::edge_histogram() {
  vector<pair<int, int>> histogram;
  edge_length_histogram(histogram);

  ostringstream oss;
  oss << "Bucket, Edges, Non-Edges" << std::endl;
  FOR(i, histogram.size()) {
    pair<int, int> bucket = histogram[i];
    oss << i << "," << bucket.first << "," << bucket.second << std::endl;
  }
  return oss.str();
}

void Embedding::edge_length_histogram(vector<std::pair<int, int>>& histogram) {
  int number_of_buckets = 100;
  // double max_distance = 2 * R;

  double max_distance = 0.0;
  FOR(i, G.n - 1) {
    for (int j = i + 1; j < G.n; j++) {
      double vertex_distance =
          Force::hyperbolic_distance_between_coordinates(points[i], points[j]);
      if (vertex_distance > max_distance) {
        max_distance = vertex_distance;
      }
    }
  }

  double bucket_size = max_distance / double(number_of_buckets);
  vector<int> edges(number_of_buckets + 1, 0);
  vector<int> non_edges(number_of_buckets + 1, 0);

  int nan_count = 0;
  FOR(i, G.n) {
    FOR(j, G.n) {
      if (i < j) {
        double vertex_distance = Force::hyperbolic_distance_between_coordinates(
            points[i], points[j]);
        if (vertex_distance == vertex_distance) {
          int bucket = int(vertex_distance / bucket_size);
          if (G.adjacent(i, j)) {
            edges[bucket] += 1;
          } else {
            non_edges[bucket] += 1;
          }
        } else {
          nan_count += 1;
        }
      }
    }
  }

  std::cout << "Ignored " << nan_count << " NaNs." << std::endl;

  FOR(i, number_of_buckets) {
    pair<int, int> bucket_content(edges[i], non_edges[i]);
    histogram.push_back(bucket_content);
  }
}

double Embedding::current_expected_degree_of_vertex(int vertex) {
  double R = this->max_radius();
  double vertex_radius = points[vertex].length();
  double first_factor = 2.0 / M_PI * exp((R - vertex_radius) / 2.0);

  double radii_sum = 0.0;
  FOR(i, G.n) {
    if (i != vertex) {
      double other_vertex_radius = points[i].length();
      radii_sum += exp(-(other_vertex_radius / 2.0));
    }
  }

  return first_factor * radii_sum;
}

Coordinate Embedding::fitting_plane() const {
  int n = G.n;
  Matrix<double, 3, Dynamic> point_cloud_matrix;

  FOR(i, n) {
    vector<double> point = {points[i].x, points[i].y, points[i].z};
    Vector3d column(point.data());
    point_cloud_matrix.conservativeResize(point_cloud_matrix.rows(),
                                          point_cloud_matrix.cols() + 1);
    point_cloud_matrix.col(point_cloud_matrix.cols() - 1) = column;
  }

  /**
   *  Perform a singular value decomposition of the point cloud matrix.
   */
  Eigen::JacobiSVD<Matrix<double, 3, Dynamic>> svd(
      point_cloud_matrix, Eigen::ComputeFullU | Eigen::ComputeFullV);
  VectorXd singular_values = svd.singularValues();
  vector<double> singulars(
      singular_values.data(),
      singular_values.data() + singular_values.rows() * singular_values.cols());

  /**
   *  Find the smallest singular value.
   */
  int index_of_minimal_singular = 0;
  double minimal_singular = std::numeric_limits<double>::max();
  FOR(i, singulars.size()) {
    if (singulars[i] < minimal_singular) {
      minimal_singular = singulars[i];
      index_of_minimal_singular = i;
    }
  }

  VectorXd minimal_left_singular_vector =
      svd.matrixU().col(index_of_minimal_singular);
  vector<double> plane_normal_vector(
      minimal_left_singular_vector.data(),
      minimal_left_singular_vector.data() +
          minimal_left_singular_vector.rows() *
              minimal_left_singular_vector.cols());

  return Coordinate(plane_normal_vector[0], plane_normal_vector[1],
                    plane_normal_vector[2]);
}

Coordinate Embedding::euclidean_centroid() const {
  Coordinate centroid = Coordinate();
  FOR(i, points.size()) { centroid += points[i]; }

  centroid /= (double)points.size();
  return centroid;
}

Coordinate Embedding::approximate_centroid() const {
  int min_degree = G.edges[G.n - 1].size();
  int first_index_of_degree_one_vertices = 0;
  for (int i = G.n - 1; i > 0; i--) {
    if (G.edges[i].size() > min_degree) {
      first_index_of_degree_one_vertices = i;
      break;
    }
  }

  /**
   *  Randomly sampling degree-1 vertices and choosing the pair with the
   *  largest distance.
   */
  std::random_device rd;   // obtain a random number from hardware
  std::mt19937 eng(rd());  // seed the generator
  std::uniform_int_distribution<> distr(first_index_of_degree_one_vertices,
                                        G.n - 1);  // define the range

  double max_distance = 0.0;
  Coordinate max_dist_point_1 = Coordinate();
  Coordinate max_dist_point_2 = Coordinate();

  for (int k = 0; k < 20; k++) {
    int vertex_1 = distr(eng);
    int vertex_2 = distr(eng);

    Coordinate point_1 = points[vertex_1];
    Coordinate point_2 = points[vertex_2];

    double distance =
        Force::hyperbolic_distance_between_coordinates(point_1, point_2);

    if (distance > max_distance) {
      max_distance = distance;
      max_dist_point_1 = point_1;
      max_dist_point_2 = point_2;
    }
  }

  if (max_distance == 0.0) {
    return Coordinate();
  }

  Coordinate max_point_1_spherical = max_dist_point_1.cartesian_to_spherical();
  Coordinate max_point_2_spherical = max_dist_point_2.cartesian_to_spherical();

  Coordinate approximate_centroid =
      HyperbolicSpace::weighted_centroid_between_coordinates(
          max_point_1_spherical, 1.0, max_point_2_spherical, 1.0);
  return approximate_centroid.spherical_to_cartesian();
}

void Embedding::center() {
  Coordinate apprximate_centroid = approximate_centroid();
  move_vector_to_origin(apprximate_centroid);
}

void Embedding::move_vector_to_origin(const Coordinate& vector) {
  if (fabs(vector.x) > 0.0000000001 || fabs(vector.y) > 0.0000000001 ||
      fabs(vector.z) > 0.0000000001) {
    FOR(i, points.size()) {
      Coordinate old_coordinate_spherical = points[i].cartesian_to_spherical();
      Coordinate new_coordinate =
          HyperbolicSpace::translate_point_by_vector(points[i], vector);
      points[i] = new_coordinate;
    }
  }
}

double Embedding::distance_from_point_to_plane_with_normal(Coordinate& point,
                                                           Coordinate& normal) {
  double denominator =
      sqrt(normal.x * normal.x + normal.y * normal.y + normal.z * normal.z);

  double enumerator =
      normal.x * point.x + normal.y * point.y + normal.z * point.z;

  return abs(enumerator) / denominator;
}

double Embedding::average_distance_of_embedding_to_plane_with_normal(
    Coordinate& normal) {
  double total_distance = 0.0;
  FOR(i, G.n) {
    total_distance +=
        distance_from_point_to_plane_with_normal(points[i], normal);
  }

  return total_distance / (double)G.n;
}

bool Embedding::is_stable() const { return is_stable_with_threshold(0.25); }

bool Embedding::is_stable_with_threshold(double threshhold) const {
  int iterations = potential_history.size();
  if (iterations > 10) {
    bool descending = true;
    for (int i = 1; i < 9; i++) {
      if (potential_history[iterations - i] >=
          potential_history[iterations - (i + 1)]) {
        descending = false;
        break;
      }
    }

    if (descending) {
      double previous_potential = potential_history[iterations - 2];
      double current_potential = potential_history[iterations - 1];

      double potential_difference = previous_potential - current_potential;
      if (potential_difference < threshhold) {
        return true;
      }
    }

    if (iterations > 20) {
      /**
       * Checking for oscillations.
       */
      double max_difference = 0.0;
      for (int i = 1; i < 19; i++) {
        double difference = fabs(potential_history[iterations - i] -
                                 potential_history[iterations - (i + 1)]);
        if (difference > max_difference) {
          max_difference = difference;
        }
      }

      if (max_difference < 1.0) {
        return true;
      }
    }
  }

  return false;
}

bool Embedding::is_radially_stable() const {
  int iterations = radial_potential_history.size();
  if (iterations >= 2) {
    double previous_potential = radial_potential_history[iterations - 2];
    double current_potential = radial_potential_history[iterations - 1];

    if (current_potential < previous_potential) {
      double potential_difference = previous_potential - current_potential;
      if (potential_difference < 0.25) {
        return true;
      }
    }
  }

  return false;
}

void Embedding::convert_to_2D_embedding() {
  Coordinate fitting_plane_normal = fitting_plane();
  fitting_plane_normal /= fitting_plane_normal.length();

  Coordinate z_axis = Coordinate(0.0, 0.0, 1.0);

  int n = G.n;
  FOR(i, n) {
    /**
     *  Now we add the force that pulls the vertex towards the plane.
     *  Therefore, we get the coordinate of the vertex on the plane.
     *
     *  Meaning we take the radius of the vertex and project its angular
     *  coordinate onto the plane.
     */
    Coordinate rotation_axis_onto_plane =
        Coordinate::normal_vector_for_plane_from_points(
            Coordinate(), fitting_plane_normal, points[i]);

    rotation_axis_onto_plane /= rotation_axis_onto_plane.length();
    Coordinate point_on_normal_vector =
        fitting_plane_normal * points[i].length();

    Coordinate point_on_plane = point_on_normal_vector.coordinate_by_rotation(
        rotation_axis_onto_plane, M_PI_2);

    /**
     *  Now all the points are on the same plane. It remains to rotated the
     *  plane such that it is equal to the xy-plane.
     */

    Coordinate rotation_axis_onto_x_y_plane =
        Coordinate::normal_vector_for_plane_from_points(Coordinate(), z_axis,
                                                        fitting_plane_normal);

    double rotation_angle =
        atan2(z_axis.cross_product(fitting_plane_normal).length(),
              z_axis.dot_product(fitting_plane_normal));

    Coordinate new_coordinate = point_on_plane.coordinate_by_rotation(
        rotation_axis_onto_x_y_plane, M_PI - rotation_angle);

    points[i] = new_coordinate;
  }
}

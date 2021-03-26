/**
 * The threege tool is used to embed a graph into the hyperbolic plane
 * by first embedding it into the 3-dimensional hyperbolic space and
 * forcing all vertices onto the plane afterwards.
 */

#include <gflags/gflags.h>
#include <glog/logging.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>

#include "NLEHelper.h"
#include "constantforce.h"
#include "coordinate.h"
#include "embedding.h"
#include "gnuplot_i.h"
#include "graph.h"
#include "surfaceforce.h"
// #include "force.hpp"
#include "hyperbolicspace.h"

#define FOR(i, n) for (int(i) = 0; (i) < (n); (i)++)

/**
 * Flags
 */

DEFINE_string(graph, "", "The input file containing the graph to embed.");
DEFINE_bool(info, true, "Print information about the graph.");
DEFINE_int32(seed, 0, "The seed used to generate random numbers.");
DEFINE_string(embed, "",
              "If set, embeds the graph and writes the resulting embedding to "
              "the specified file.");
DEFINE_string(print, "", "Writes the graph to the specified file.");
DEFINE_string(draw, "",
              "Saves the embedding as SVG file to the specified path.");
DEFINE_string(embedding, "", "Reads an embedding from the passed path.");
DEFINE_string(compare, "",
              "Compare the calculated embedding to the one saved in the passed "
              "file path. (phi-phi / r-r)");
DEFINE_string(plot, "",
              "Saves the plots representing the quality of the embedding to "
              "the passed filepath.");

DEFINE_int32(degree, 5,
             "Only vertices whose degree is larger than this value will have "
             "node label in the drawn SVG file.");

DEFINE_string(names, "",
              "Path to a file that maps the node indices to strings. When "
              "drawing the graph as SVG these names will be used instead of "
              "the indices.");

DEFINE_string(edgeHisto, "",
              "File where the edge length histogram should be plotted to.");

DEFINE_string(
    phase1Out, "",
    "Where the embedding after the first phase should be written to.");

// Flags that represent parameters to tweak the embedding process.

DEFINE_double(stable, 0.1,
              "Determines the threshold for when an embedding counts as "
              "stable or not.");

DEFINE_double(repulsiveForceCoefficient, 0.03,
              "Determines how repulsive forces should to be scaled. Usually, "
              "the should be much weaker than the attractive forces, so a "
              "small value of 0.03 was set as default.");

DEFINE_double(
    repulsiveForceCoefficientSecondPhase, 0.03,
    "Determines how repulsive forces should to be scaled in the "
    "second phase. Usually, the should be much weaker than the "
    "attractive forces, so a small value of 0.03 was set as default.");

DEFINE_double(attractiveExponent, 1.0,
              "Determines the exponent of the attractive force.");

DEFINE_double(repulsiveExponent, 0.5,
              "Determines the exponent of the repulsive force.");

DEFINE_double(velocity, 0.5,
              "The velocity of the vertices which determines how much of the "
              "previous forces are again part of the current force.");

DEFINE_double(initialTemperature, 0.5,
              "The initial temperature of the system.");

DEFINE_double(temperatureDecreaseFactor, 0.975,
              "The factor that the temperature will be multiplied with, after "
              "an iteration. E.g. a factor of 0.5 will halve the temperature "
              "in each iteration.");


void print_edge_stats_of_embedding(const Embedding& embedding) {
  std::cout << "Determining edge stats..."
            << "\n";
  double number_of_edges = 0.0;
  double number_of_non_edges = 0.0;
  double long_edges = 0.0;
  double short_non_edges = 0.0;
  double R = embedding.max_radius();
  for (int i = 0; i < embedding.G.n - 1; i++) {
    Coordinate i_coord = embedding.points[i];
    for (int j = i + 1; j < embedding.G.n; j++) {
      Coordinate j_coord = embedding.points[j];
      double edge_length =
          Force::hyperbolic_distance_between_coordinates(i_coord, j_coord);

      double difference_to_desired = R - edge_length;

      if (embedding.G.adjacent(i, j)) {
        if (difference_to_desired < 0.0) {
          long_edges += difference_to_desired * difference_to_desired;
        }
        number_of_edges += 1.0;
      } else {
        if (difference_to_desired > 0.0) {
          short_non_edges += difference_to_desired * difference_to_desired;
        }
        number_of_non_edges += 1.0;
      }
    }
  }

  std::cout << "Edges too long: " << long_edges / number_of_edges << "\n";
  std::cout << "Non-Edges too short: " << short_non_edges / number_of_non_edges
            << "\n";
}

void plot_edge_histogram_of_embedding_to_file(Embedding& embedding,
                                              string file) {
  std::cout << "Plotting edge length histogram to: " << file << "\n";
  vector<std::pair<int, int>> histogram;
  embedding.edge_length_histogram(histogram);

  vector<int> indices(histogram.size());
  vector<double> edge_buckets(histogram.size());
  vector<double> non_edge_buckets(histogram.size());

  double number_of_edges = 0.0;
  double number_of_non_edges = 0.0;

  for (int i = 0; i < (int)histogram.size(); i++) {
    std::pair<int, int> bucket = histogram[i];

    double number_of_edges_in_bucket = (double)bucket.first;
    number_of_edges += number_of_edges_in_bucket;
    edge_buckets[i] = number_of_edges_in_bucket;

    double number_of_non_edges_in_bucket = (double)bucket.second;
    number_of_non_edges += number_of_non_edges_in_bucket;
    non_edge_buckets[i] = number_of_non_edges_in_bucket;

    indices[i] = i;
  }

  for (int i = 0; i < (int)histogram.size(); i++) {
    edge_buckets[i] /= number_of_edges;
    non_edge_buckets[i] /= number_of_non_edges;
  }

  Gnuplot gp("Edge Length Histogram");
  gp.savetops(file);

  // set range
  gp.set_xrange(0, indices.size());
  // gp.set_yrange(0, 1.0);

  gp.set_style("lines lw 5 lc rgb \"green\"")
      .plot_xy(indices, edge_buckets, "Edges");

  gp.set_style("lines lw 5 lc rgb \"red\"")
      .plot_xy(indices, non_edge_buckets, "Non Edges");

  /**
   * Now we determine the error we made using this embedding.
   *
   * This works as follows. For each bucket we check which of
   * the two (edges or non-edges) is smaller and add the bucket
   * value to the corresponding error term.
   */

  double edges_error = 0.0;
  double non_edges_error = 0.0;

  for (int i = 0; i < (int)histogram.size(); i++) {
    if (edge_buckets[i] > non_edge_buckets[i]) {
      non_edges_error += non_edge_buckets[i];
    } else {
      edges_error += edge_buckets[i];
    }
  }

  std::cout << "Histogram Error Edges: " << edges_error << "\n";
  std::cout << "Histogram Error Non Edges: " << non_edges_error << "\n";
}

/**
 * Main procedure
 */
int main(int argc, char* argv[]) {
  /**
   * Evaluating the flags.
   */
  google::ParseCommandLineFlags(&argc, &argv, true);
  google::InitGoogleLogging(argv[0]);

  FLAGS_log_dir = "/tmp/";

  std::ios_base::sync_with_stdio(false);

  std::srand(FLAGS_seed);

  if (FLAGS_graph.empty()) {
    LOG(WARNING) << "No input graph specified!";
    return EXIT_FAILURE;
  }

  unordered_map<std::string, int> node_labels;
  vector<int> new_permutation;
  auto G = Graph::fromFile(FLAGS_graph, &node_labels);
  G.sortByDegrees(&new_permutation);

  for (auto label_pair : node_labels) {
    string label = label_pair.first;
    int old_index = label_pair.second;

    node_labels[label] = new_permutation[old_index];
  }

  if (FLAGS_info) {
    std::cout << "# Graph: " << FLAGS_graph << "\n";
    std::cout << "# n = " << G.n << "\n";
    std::cout << "# m = " << G.numEdges() << "\n";
    std::cout << "# avg. deg. = " << G.averageDegree() << "\n";
    std::cout << "# beta = " << G.powerLawExponent() << "\n";
    std::cout << "# clustering = " << G.clusterCoeff() << "\n";

    double T, n_orig, m_orig, alpha, R;
    NLEHelper::estimateHyperbolicParameters(G, &T, &n_orig, &m_orig, &alpha, &R,
                                            nullptr);

    std::cout << "# Hyperbolic Parameters: n = " << n_orig << ", m = " << m_orig
              << ", R = " << R << ", T = " << T << ", alpha = " << alpha
              << "\n";
  }

  int n = G.n;
  int total_iterations = 0;
  double T, n_orig, m_orig, alpha, R;
  vector<double> estimated_radii(n);
  NLEHelper::estimateHyperbolicParameters(G, &T, &n_orig, &m_orig, &alpha, &R,
                                          &estimated_radii);

  if (FLAGS_print.length() > 0) {
    std::cout << "Writing graph to file: " << FLAGS_print << "\n";
    G.printToFile(FLAGS_print.c_str(), nullptr, false);
  }

  if (FLAGS_embed.length() > 0) {
    std::chrono::high_resolution_clock::time_point embedding_start =
        std::chrono::high_resolution_clock::now();
    double max_radius = 0.0;
    vector<Coordinate> points(n);
    FOR(i, n) {
      int seed = rand();
      Coordinate rand_coord =
          Coordinate::random_point_on_surface_of_unit_ball(seed);
      std::uniform_real_distribution<> distribution(0, R);
      int seed2 = rand();
      std::default_random_engine generator(seed2);
      double radius = distribution(generator);

      if (radius > max_radius) {
        max_radius = radius;
      }

      rand_coord *= radius;
      points[i] = rand_coord;
    }

    Embedding embedding = Embedding(G, points);

    embedding.velocity = FLAGS_velocity;

    /**
     * Before we start we determine the edge length histogram. This
     * will be overwritten later but at least we will have the error
     * terms in the console output.
     */
    if (FLAGS_edgeHisto.length() > 0) {
      plot_edge_histogram_of_embedding_to_file(embedding, FLAGS_edgeHisto);
    }

    /**
     * Phase 1: General Forces...
     */
    std::cout << "Phase 1: Radial + General forces..."
              << "\n";

    /**
     * This force is responsible for moving the nodes on the sphere surface.
     */
    SurfaceForce force = SurfaceForce(R, M_PI_4 / 2.0);
    force.temperature = FLAGS_initialTemperature;
    force.temperature_decrease = FLAGS_temperatureDecreaseFactor;
    force.repulsive_force_coefficient = FLAGS_repulsiveForceCoefficient;
    force.attractive_exponent = FLAGS_attractiveExponent;
    force.repulsive_exponent = FLAGS_repulsiveExponent;

    /**
     * This force is responsible for determining the radii.
     */
    SurfaceForce radial_force = SurfaceForce(max_radius, M_PI_4 / 2.0);
    radial_force.temperature = FLAGS_initialTemperature;
    radial_force.temperature_decrease = FLAGS_temperatureDecreaseFactor;
    radial_force.repulsive_force_coefficient = FLAGS_repulsiveForceCoefficient;

    /**
     * Start with 10 iterations of radial forces.
     */
    int initial_radial_iterations = 10;
    std::cout << initial_radial_iterations << " iterations of radial forces..."
              << "\n";

    FOR(i, initial_radial_iterations) {
      embedding.apply_radial_force_to_vertices(radial_force);

      /**
       * After applying the radial forces, we determine the new
       * maximum radius, which will be our new disk radius.
       */
      double new_max_radius = embedding.max_radius();
      force.disk_radius = new_max_radius;
    }

    std::cout << initial_radial_iterations
              << " iterations of radial forces... Done."
              << "\n";

    int current_iteration = 0;
    int max_iteration = 500;

    int radial_forces_interval = 20;
    int radial_forces_intermediate_iterations = 10;

    std::cout << "Performing at most " << max_iteration
              << " iterations of general forces and applying "
              << radial_forces_intermediate_iterations
              << " iterations of radial forces every " << radial_forces_interval
              << " iterations."
              << "\n";

    bool embedding_is_stable = false;
    // double average_percent_of_processed_cells = 0.0;
    while (!embedding_is_stable && current_iteration < max_iteration) {
      /**
       * Every "radial_forces_interval" iterations we apply radial forces.
       */
      if (current_iteration % radial_forces_interval == 0) {
        FOR(i, radial_forces_intermediate_iterations) {
          embedding.apply_radial_force_to_vertices(radial_force);

          /**
           * After applying the radial forces, we determine the new
           * maximum radius, which will be our new disk radius.
           */
          double new_max_radius = embedding.max_radius();
          force.disk_radius = new_max_radius;
        }
      }

      current_iteration++;
      total_iterations++;

      embedding.apply_force_to_vertices(force);

      /**
       * Check if the embedding is stable now.
       */
      embedding_is_stable = embedding.is_stable_with_threshold(FLAGS_stable);

      if (force.temperature * force.temperature_decrease > 0.015) {
        force.temperature *= force.temperature_decrease;
      } else {
        force.temperature = 0.015;
      }

      /**
       * We also decrease the temperature of the radial_force.
       */
      if (radial_force.temperature * radial_force.temperature_decrease >
          0.015) {
        radial_force.temperature *= radial_force.temperature_decrease;
      } else {
        radial_force.temperature = 0.015;
      }
    }

    /**
     * Stop iterating if the maximum number of iterations is reached.
     */
    if (current_iteration < max_iteration) {
      std::cout << "Stopped early (" << current_iteration << " of "
                << max_iteration
                << " iterations) since the embedding was already stable."
                << "\n";
    }

    /**
     * Save the 3D embedding after the first phase.
     */
    if (!FLAGS_phase1Out.empty()) {
      std::cout << "Saving embedding after first phase. File: "
                << FLAGS_phase1Out << "\n";
      std::ofstream out_file;
      out_file.open(FLAGS_phase1Out);
      out_file << embedding.to_string(true);
      out_file.close();
    }

    /**
     * Phase 1 End
     */

    /**
     * Phase 2: Plane Forces...
     */
    std::cout << "Phase 2: Plane forces..."
              << "\n";

    /**
     * Determine the plane that the vertices should be pulled
     * towards.
     */
    Coordinate plane_normal = embedding.fitting_plane();

    /**
     * The force that pulls the vertices towards the plane.
     */
    ConstantForce plane_force = ConstantForce(max_radius, M_PI / 3.0);
    plane_force.temperature = force.temperature;

    /**
     * The force that we used in the first phase will still by
     * applied in the second phase. However now we make the repulsive
     * forces larger, as it was too weak (on purpose) in the first
     * phase.
     */
    force.repulsive_force_coefficient =
        FLAGS_repulsiveForceCoefficientSecondPhase;

    current_iteration = 0;
    max_iteration = 500;

    std::cout << "Performing at most " << max_iteration << " iterations."
              << "\n";

    /**
     * We stop when the average distances of the vertices to the
     * plane is small enough.
     */
    double average_distance_to_plane =
        embedding.average_distance_of_embedding_to_plane_with_normal(
            plane_normal);

    while (average_distance_to_plane > 0.5 &&
           current_iteration < max_iteration) {
      current_iteration++;
      total_iterations++;

      // embedding
      //     .apply_approximate_force_to_vertices_and_flatten_towards_plane_with_normal(
      //         force, plane_force, plane_normal, FLAGS_pooling * R);
      embedding.apply_force_to_vertices_and_flatten_towards_plane_with_normal(
          force, plane_force, plane_normal);

      average_distance_to_plane =
          embedding.average_distance_of_embedding_to_plane_with_normal(
              plane_normal);

      /**
       * Update the temperature of the force
       */
      if (force.temperature * force.temperature_decrease > 0.015) {
        force.temperature *= force.temperature_decrease;
      } else {
        force.temperature = 0.015;
      }

      plane_force.temperature = force.temperature;
    }

    /**
     * Stop if we have already reached the maximum number of iterations.
     */
    if (current_iteration < max_iteration) {
      std::cout << "Stopped early (" << current_iteration << " of "
                << max_iteration
                << " iterations) since the average distance to the plane "
                   "was already small enough."
                << "\n";
    }

    /**
     * Phase 2 End
     */

    /**
     * Phase 3: Reduce to two dimensions.
     */
    std::cout << "Phase 3: Reduce to two dimensions."
              << "\n";
    embedding.convert_to_2D_embedding();

    /**
     * Phase 3 End
     */

    std::chrono::high_resolution_clock::time_point embedding_end =
        std::chrono::high_resolution_clock::now();
    auto embedding_duration =
        std::chrono::duration_cast<std::chrono::milliseconds>(embedding_end -
                                                              embedding_start)
            .count();
    std::cout << "The embedding process took: " << embedding_duration / 1000.0
              << " seconds."
              << "\n";

    std::cout << "Writing embedding to file: " << FLAGS_embed << "\n";
    std::ofstream file;
    file.open(FLAGS_embed);
    file << embedding.to_string_2D(true) << "\n";
    file.close();

    if (FLAGS_draw.length() > 0) {
      std::cout << "Drawing embedding to file: " << FLAGS_draw << "\n";
      vector<int> node_color(embedding.G.n, 1);

      if (FLAGS_names.length() > 0) {
        std::cout << "Reading node names from file: " << FLAGS_names << "\n";
        std::ifstream file;
        file.open(FLAGS_names.c_str());
        CHECK(file.good()) << "Maybe the file " << FLAGS_names
                           << " doesn't exist?";

        for (std::string line; getline(file, line);) {
          std::vector<string> string_components;
          std::stringstream ss(line);
          std::string component;

          while (ss >> component) {
            string_components.push_back(component);
            if (ss.peek() == ' ') {
              ss.ignore();
            }
          }

          string node_label = string_components[0];
          string node_name = "";

          for (int i = 1; i < string_components.size(); i++) {
            node_name += " " + string_components[i];
          }

          int node_index = node_labels[node_label];

          embedding.G.labels[node_index] = node_name;
        }
      }

      embedding.draw(FLAGS_draw.c_str(), node_color, FLAGS_degree, "", 1000,
                     1000);

      std::cout << "Embedding Log-Likelihood: " << embedding.log_likelihood()
                << "\n";
    }

    if (FLAGS_compare.length() > 0 && FLAGS_plot.length() > 0) {
      std::cout << "Comparing with embedding: " << FLAGS_compare << "\n";
      try {
        // create a new gnuplot file
        Gnuplot gp("Embedding Quality");
        gp.savetops(FLAGS_plot);

        // set range
        gp.set_xrange(0.0, 2.0 * M_PI);
        gp.set_yrange(0.0, 2.0 * M_PI);

        std::cout << "Reading embedding from: " << FLAGS_compare << "\n";
        vector<Coordinate> other_points(embedding.G.n, Coordinate());
        Embedding other_embedding = Embedding(embedding.G, other_points);
        other_embedding.read_embedding_from_file(FLAGS_compare.c_str(),
                                                 node_labels);

        std::cout << "Other Embedding Log-Likelihood: "
                  << other_embedding.log_likelihood() << "\n";

        std::cout << "Plotting data to: " << FLAGS_plot << "\n";
        vector<double> other_phi_values(embedding.G.n);
        other_embedding.get_phi_values(other_phi_values);

        vector<double> this_phi_values(embedding.G.n);
        embedding.get_phi_values(this_phi_values);

        gp.set_style("points pointtype 5 ps 0.75 lc rgb \"blue\"")
            .plot_xy(other_phi_values, this_phi_values, "Phi-Phi Plot");

        // start a new picture
        gp.reset_plot();

        double R = embedding.max_radius();

        gp.set_xrange(0.0, 1.5 * R);
        gp.set_yrange(0.0, 1.5 * R);

        vector<double> other_r_values(embedding.G.n);
        other_embedding.get_r_values(other_r_values);

        vector<double> this_r_values(embedding.G.n);
        embedding.get_r_values(this_r_values);

        gp.set_style("points pointtype 5 ps 0.75 lc rgb \"blue\"")
            .plot_xy(other_r_values, this_r_values, "R-R Plot");

      } catch (GnuplotException ge) {
        cout << ge.what() << endl;
      }
    }

    if (FLAGS_edgeHisto.length() > 0) {
      plot_edge_histogram_of_embedding_to_file(embedding, FLAGS_edgeHisto);
    }
  }

  if (FLAGS_embedding.length() > 0) {
    std::cout << "Reading embedding from file: " << FLAGS_embedding << "\n";
    vector<Coordinate> points(G.n, Coordinate());
    Embedding embedding = Embedding(G, points);
    embedding.read_embedding_from_file(FLAGS_embedding.c_str(), node_labels);

    if (FLAGS_draw.length() > 0) {
      std::cout << "Drawing embedding to file: " << FLAGS_draw << "\n";
      vector<int> node_color(embedding.G.n, 1);
      embedding.draw(FLAGS_draw.c_str(), node_color, FLAGS_degree, "", 1000,
                     1000);

      // if (FLAGS_show) {
      //   string open_drawing_command = "open " + FLAGS_draw;
      //   system(open_drawing_command.c_str());
      // }
    }

    if (FLAGS_compare.length() > 0 && FLAGS_plot.length() > 0) {
      std::cout << "Comparing with embedding: " << FLAGS_compare << "\n";
      try {
        // create a new gnuplot file
        Gnuplot gp("Embedding Quality");
        gp.savetops(FLAGS_plot);

        // set range
        gp.set_xrange(0.0, 2.0 * M_PI);
        gp.set_yrange(0.0, 2.0 * M_PI);

        std::cout << "Reading embedding from: " << FLAGS_compare << "\n";
        vector<Coordinate> other_points(embedding.G.n, Coordinate());
        Embedding other_embedding = Embedding(embedding.G, other_points);
        other_embedding.read_embedding_from_file(FLAGS_compare.c_str(),
                                                 node_labels);

        std::cout << "Other Embedding Log-Likelihood: "
                  << other_embedding.log_likelihood() << "\n";

        // if (FLAGS_edgeStats) {
        //   print_edge_stats_of_embedding(other_embedding);
        // }

        std::cout << "Plotting data to: " << FLAGS_plot << "\n";
        vector<double> other_phi_values(embedding.G.n);
        other_embedding.get_phi_values(other_phi_values);

        vector<double> this_phi_values(embedding.G.n);
        embedding.get_phi_values(this_phi_values);

        gp.set_style("points pointtype 5 ps 0.75 lc rgb \"blue\"")
            .plot_xy(other_phi_values, this_phi_values, "Phi-Phi Plot");

        // start a new picture
        gp.reset_plot();

        double R = embedding.max_radius();

        gp.set_xrange(0.0, 1.5 * R);
        gp.set_yrange(0.0, 1.5 * R);

        vector<double> other_r_values(embedding.G.n);
        other_embedding.get_r_values(other_r_values);

        vector<double> this_r_values(embedding.G.n);
        embedding.get_r_values(this_r_values);

        gp.set_style("points pointtype 5 ps 0.75 lc rgb \"blue\"")
            .plot_xy(other_r_values, this_r_values, "R-R Plot");

        // if (FLAGS_show) {
        //   string open_plot_command = "open " + FLAGS_plot;
        //   system(open_plot_command.c_str());
        // }

      } catch (GnuplotException ge) {
        cout << ge.what() << endl;
      }
    }

    if (FLAGS_edgeHisto.length() > 0) {
      plot_edge_histogram_of_embedding_to_file(embedding, FLAGS_edgeHisto);
    }
  }

  return 0;
}

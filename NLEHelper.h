//
//  NLEHelper.hpp
//  graphlib
//
//  Created by Anton Krohmer on 11.11.15.
//
//

#pragma once

#include <fstream>
#include <functional>
#include <iostream>

#include "graph.h"

typedef double (*Integrable)(double* x_array, size_t dim, void* params);

class NLEHelper {
 public:
  static void estimateHyperbolicParameters(const Graph& H, double* T, double* n,
                                           double* m, double* alpha, double* R,
                                           vector<double>* radial_coords);

  static void estimateHyperbolicParameters(const Graph& H, std::ostream& os);

  static void estimateHyperbolicParameters(std::istream& is, double* T,
                                           double* n, double* m, double* alpha,
                                           double* R,
                                           vector<double>* radial_coords);

  static double integrate(Integrable f, int dim, double upper_bds[],
                          double lower_bds[], void* params);

  // Estimates the original N.
  double computeFitnessN(double n_orig);
  // Estimates the parameter R.
  double computeFitnessR(double R);

  // Computes the difference between the estimated number of edges of a node
  // at estimated_coord and num_edges.
  double radialCoordFitness(int num_edges, double estimated_coord);

 private:
  NLEHelper(const Graph& H, double T) : H(H), T(T) {
    m = 0;
    for (const auto& node : H.edges) m += node.size();
    m /= 2;
  }

  // Uses the Newman method to retrieve the power-law.
  void estimatePowerLaw();

  // Uses a CMAES to estimate the graph parameters.
  void estimateGlobalParameters();

  void estimateRadialCoordinates();

  // Optimizes the given function by binary search on to_opt until a
  // fixed optimum.
  // Left and right are bounds on where the optimum is.
  void optimize(double* to_opt, const std::function<double(double)>& fitness,
                double left = 0, double right = 0, double precision = 0.01);

  const Graph& H;
  double T;
  // Number of edges
  int m;

  // Estimated parameters
  double n_orig, m_orig;
  double alpha, R;

  vector<double>* radial_coords;
};

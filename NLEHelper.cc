//
//  NLEHelper.cpp
//  graphlib
//
//  Created by Anton Krohmer on 11.11.15.
//
//

#include "NLEHelper.h"

#include <glog/logging.h>
#include <gsl/gsl_monte_vegas.h>
#include <gsl/gsl_sf_gamma.h>

#include <limits>

#include "powerlawCommon.h"
#include "random.h"

#define FOR(i, n) for (int(i) = 0; (i) < (n); (i)++)
#define FORB(i, a, n) for (int(i) = (a); (i) < (n); (i)++)
#define FOREACH(it, c) \
  for (__typeof((c).begin()) it = (c).begin(); it != (c).end(); ++it)
#define PB emplace_back
#define MP make_pair

using std::placeholders::_1;
using std::placeholders::_2;
using std::placeholders::_3;

double NLEHelper::integrate(Integrable f, int dim, double upper_bds[],
                            double lower_bds[], void* params) {
  double res, err;
  gsl_monte_function F;
  F.f = f;
  F.params = params;
  F.dim = dim;

  gsl_monte_vegas_state* s = gsl_monte_vegas_alloc(dim);

  gsl_monte_vegas_integrate(&F, lower_bds, upper_bds, dim, 2000, rng(), s, &res,
                            &err);

  while (fabs(gsl_monte_vegas_chisq(s) - 1.0) > 0.5 && err > 1e-80) {
    //        LOG(INFO) << "Main iteration!";
    gsl_monte_vegas_integrate(&F, lower_bds, upper_bds, dim, 20000, rng(), s,
                              &res, &err);
    //            printf("result = % .6f sigma = % .6f "
    //                   "chisq/dof = %.1f\n", res, err, gsl_monte_vegas_chisq
    //                   (s));
    CHECK(!std::isnan(res));
  }

  gsl_monte_vegas_free(s);

  // LOG(INFO) << "res=" << res << ", err=" << err;
  return res;
}

namespace {
struct exp_deg_params {
  double alpha;
  int n;
  double R;
  double T;
  double r_v;
};

double exp_deg(double x[], size_t dim, void* p) {
  struct exp_deg_params* fp = (struct exp_deg_params*)p;

  CHECK_EQ(dim, 2);
  double r_1 = x[0];
  double phi = x[1];

  // prevent numerical errors
  double dist =
      cosh(r_1) * cosh(fp->r_v) - sinh(r_1) * sinh(fp->r_v) * cos(phi);
  if (dist < 1) dist = 1;

  double ret = (1 / M_PI) * (fp->n - 1) * fp->alpha * sinh(fp->alpha * r_1) /
               ((cosh(fp->alpha * fp->R) - 1) *
                (1 + exp((1 / (2 * fp->T)) * (acosh(dist) - fp->R))));

  return ret;
}

struct avg_deg_params {
  double alpha;
  double n;
  double R;
  double T;  // n is double for numerical
             // solving
};

double avg_deg(double x[], size_t dim, void* p) {
  struct avg_deg_params* fp = (struct avg_deg_params*)p;

  CHECK_EQ(dim, 3);
  double r_1 = x[0];
  double r_2 = x[1];
  double phi = x[2];

  // prevent numerical errors
  double dist = cosh(r_1) * cosh(r_2) - sinh(r_1) * sinh(r_2) * cos(phi);
  if (dist < 1) dist = 1;

  double ret = (1 / M_PI) * (fp->n - 1) * fp->alpha * fp->alpha *
               sinh(fp->alpha * r_1) * sinh(fp->alpha * r_2) /
               ((cosh(fp->alpha * fp->R) - 1) * (cosh(fp->alpha * fp->R) - 1) *
                (1 + exp((1 / (2 * fp->T)) * (acosh(dist) - fp->R))));

  return ret;
}

double small_comp(double x[], size_t dim, void* p) {
  struct avg_deg_params* fp = (struct avg_deg_params*)p;

  CHECK_EQ(dim, 1);
  double r = x[0];

  double p1 =
      exp(-fp->n * exp(-(fp->alpha - 0.5) * fp->R - (1 - fp->alpha) * r) * 2 *
          fp->alpha / (M_PI * (fp->alpha - 0.5)));
  double p2 =
      1 +
      (fp->n * 2 * fp->alpha *
       (exp(-0.5 * r) - exp(-(fp->alpha - 0.5) * fp->R - (1 - fp->alpha) * r)) /
       (M_PI * (fp->alpha - 0.5)));
  double p3 = fp->alpha * sinh(fp->alpha * r) / (cosh(fp->alpha * fp->R) - 1);
  double ret = fp->n * p1 * p2 * p3;

  return ret;
}

//  double small_comp_edges(double x[], size_t dim, void * p) {
//    struct avg_deg_params * fp = (struct avg_deg_params *)p;
//
//    CHECK_EQ(dim, 1);
//    double r = x[0];
//
//    double p1 = exp(-fp->n * exp(-(fp->alpha - 0.5) * fp->R - (1-fp->alpha) *
//    r)
//                    * 2 * fp->alpha / (M_PI * (fp->alpha - 0.5)));
//    double p2 = (fp->n * 2 * fp->alpha * (exp(-0.5*r)
//         - exp(-(fp->alpha - 0.5) * fp->R - (1-fp->alpha) * r)
//         ) / (M_PI * (fp->alpha - 0.5)));
//    double p3 = fp->alpha * sinh(fp->alpha * r) / (cosh(fp->alpha * fp->R) -
//    1); double ret = fp->n * p1 * p2 * p3;
//
//    return ret;
//  }
}  // namespace

double NLEHelper::computeFitnessR(double R) {
  // works almost perfectly for exact alpha
  struct avg_deg_params p = {alpha, .0 + n_orig, R, T};

  double xl2[] = {0, 0, 0};
  double xu2[] = {R, R, M_PI};

  double exp_avg_deg = integrate(&avg_deg, 3, xu2, xl2, &p);

  //  double xl3[] = {R/2};
  //  double xu3[] = {R};
  //  double missing_edges = integrate(&small_comp_edges, 1, xu3, xl3, &p);
  //
  //  LOG(INFO) << "Missing edges=" << missing_edges;

  double fitness2 =
      (exp_avg_deg * n_orig - (2 * m)) * (exp_avg_deg * n_orig - (2 * m));
  //
  //  LOG(INFO) << "R=" << R << ", avgdeg=" << exp_avg_deg << ", fit=" <<
  //  fitness2;

  return fitness2;
}

double NLEHelper::computeFitnessN(double n_orig) {
  if (n_orig < H.n) return std::numeric_limits<double>::max();

  struct avg_deg_params p = {alpha, .0 + n_orig, R, T};

  //  Hyperbolic* H2 = HyperbolicLinear::linearSampling(n_orig, R, alpha, H.T);
  //  Hyperbolic* H2_giant = H2->giantSubgraph();
  //  int missing_nodes = n_orig - H2_giant->n;

  double xl[] = {R / 2};
  double xu[] = {R};
  int missing_nodes = (int)integrate(&small_comp, 1, xu, xl, &p);

  double fitness1 =
      (n_orig - missing_nodes - H.n) * (n_orig - missing_nodes - H.n);

  //  LOG(INFO) << "n_orig=" << n_orig << ", missing=" << missing_nodes
  //            << ", fitness=" << fitness1;
  //
  //  delete H2;
  //  delete H2_giant;

  return fitness1;
}

void NLEHelper::optimize(double* to_opt,
                         const std::function<double(double)>& fitness,
                         double left, double right, double precision) {
  // begin by finger search
  while (fitness(*to_opt / 2) < fitness(*to_opt)) *to_opt /= 2;
  while (fitness(*to_opt * 2) < fitness(*to_opt)) *to_opt *= 2;

  // now optimum is between to_opt/2 and 2to_opt
  left = std::max(*to_opt / 2, left);
  if (right > 0)
    right = std::min(*to_opt * 2, right);
  else
    right = *to_opt * 2;

  // use golden section search
  double golden = (sqrt(5) - 1) / 2;

  double x1 = right - golden * (right - left);
  double x2 = left + golden * (right - left);
  while (std::abs(left - right) > precision * (left + right)) {
    double fitx1 = fitness(x1);
    double fitx2 = fitness(x2);
    if (fitx1 < fitx2) {
      right = x2;
      x2 = x1;
      x1 = right - golden * (right - left);
    } else {
      left = x1;
      x1 = x2;
      x2 = left + golden * (right - left);
    }
  }

  *to_opt = (right + left) / 2;
}

void NLEHelper::estimateGlobalParameters() {
  // Use a binary search approach for indivudal variable, combine with
  // iterative improvement of other variables.

  vector<int> hist = H.degHisto();
  int missing_zeros = hist[1] - (hist[2] - hist[1]);

  n_orig = H.n + (missing_zeros > 0 ? missing_zeros : 0);
  m_orig = m;
  //  R = 2*log(H.n);

  //  optimize(&n_orig,
  //           std::bind(&NLEHelper::computeFitnessN, this, _1), H.n, 0,
  //           0.0001);
  //  optimize(&R,
  //           std::bind(&NLEHelper::computeFitnessR, this, _1), 0, 0, 0.0001);
  //  LOG(INFO) << "Integral R=" << R;

  R = 2 * log(8 * n_orig * alpha * alpha * T /
              (sin(M_PI * T) * (2. * m_orig / n_orig) * (2 * alpha - 1) *
               (2 * alpha - 1)));

  //  LOG(INFO) << "Analystical R=" << R;
}

void NLEHelper::estimatePowerLaw() {
  plfit::VectorType degs;
  FOR(i, H.n)
  degs.PB(H.edges[i].size());
  plfit::VectorType results;
  plfit::Powerlaw::SingleFit(degs, results, false, false, 1.5, 0.01, 4);

  alpha = (results[0] - 1) / 2;
  if (alpha <= 0.55) {
    LOG(WARNING) << "alpha estimated at " << alpha
                 << ", too low for embedding. Using alpha=0.55 instead.";
    alpha = 0.55;
  } else if (alpha >= 0.95) {
    LOG(WARNING) << "alpha estimated at " << alpha
                 << ", too high for embedding. Using alpha=0.95 instead.";
    alpha = 0.95;
  }

  //  LOG(INFO) << "estimated alpha=" << alpha;
  //  LOG(INFO) << "Xmin=" << results[1];
  //  LOG(INFO) << "Log-likelihood=" << results[2];
  //  LOG(INFO) << "Karl's alpha=" << (H.powerLawExponent()-1)/2;
}

void NLEHelper::estimateHyperbolicParameters(const Graph& G, std::ostream& os) {
  double R, T, alpha;
  double n_orig, m_orig;
  vector<double> est_r(G.n);
  estimateHyperbolicParameters(G, &T, &n_orig, &m_orig, &alpha, &R, &est_r);

  os << T << std::endl;
  os << n_orig << std::endl;
  os << m_orig << std::endl;
  os << alpha << std::endl;
  os << R << std::endl;

  for (int i = 0; i < G.n; ++i) {
    os << i << " " << est_r[i] << std::endl;
  }
}

void NLEHelper::estimateHyperbolicParameters(std::istream& is, double* T,
                                             double* n, double* m,
                                             double* alpha, double* R,
                                             vector<double>* radial_coords) {
  is >> *T;
  is >> *n;
  is >> *m;
  is >> *alpha;
  is >> *R;

  LOG(INFO) << "Estimated data:";
  LOG(INFO) << "alpha  = " << *alpha;
  LOG(INFO) << "T      = " << *T;
  LOG(INFO) << "n_orig = " << *n;
  LOG(INFO) << "m_orig = " << *m;
  LOG(INFO) << "R      = " << *R;

  int v;
  while (is >> v) {
    double r;
    is >> r;
    (*radial_coords)[v] = r;
  }
}

void NLEHelper::estimateHyperbolicParameters(const Graph& G, double* T,
                                             double* n, double* m,
                                             double* alpha, double* R,
                                             vector<double>* radial_coords) {
  // Setting T to a small arbitrary constant seems to produce good enough
  // results
  *T = 0.1;
  NLEHelper nle(G, *T);
  nle.radial_coords = radial_coords;

  nle.estimatePowerLaw();
  nle.estimateGlobalParameters();

  // TODO fix
  //  nle.n_orig = 18000;
  //  nle.m_orig = 175927;
  //  nle.alpha = 0.51;
  //  nle.R = 25;

  //  const Hyperbolic& H = dynamic_cast<const Hyperbolic &>(G);
  //  nle.n_orig = H.n;
  //  nle.m_orig = H.averageDegree() * H.n / 2;
  //  nle.alpha = H.alpha;
  //  nle.R = H.R;

  *n = nle.n_orig;
  *m = nle.m_orig;
  *alpha = nle.alpha;
  *R = nle.R;

  LOG(INFO) << "Estimated data:";
  LOG(INFO) << "alpha  = " << *alpha;
  LOG(INFO) << "T      = " << *T;
  LOG(INFO) << "n_orig = " << *n;
  LOG(INFO) << "m_orig = " << *m;
  LOG(INFO) << "R      = " << *R;

  //  FOR(i,H.n)
  //    radial_coords->at(i) = H.pts[i].r;
  if (radial_coords != nullptr) nle.estimateRadialCoordinates();
}

// void NLEHelper::estimateHyperbolicParameters(const Graph& H, double T,
//                                              double *n, double *m,
//                                              double *alpha, double *R) {
//   NLEHelper nle(H, T);
//   // nle.radial_coords = radial_coords;

//   nle.estimatePowerLaw();
//   nle.estimateGlobalParameters();

// //  nle.n_orig = 35685;
// //  nle.m_orig = 175927;
// //  nle.alpha = 0.55;
// //  nle.R = 27;

//   *n = nle.n_orig;
//   *m = nle.m_orig;
//   *alpha = nle.alpha;
//   *R = nle.R;

//   LOG(INFO) << "Estimated data:";
//   LOG(INFO) << "alpha  = " << *alpha;
//   LOG(INFO) << "n_orig = " << *n;
//   LOG(INFO) << "m_orig = " << *m;
//   LOG(INFO) << "R      = " << *R;

// //  nle.estimateRadialCoordinates();
// }

double NLEHelper::radialCoordFitness(int num_edges, double estimated_coord) {
  double xl[] = {0, 0};
  double xu[] = {R, M_PI};
  struct exp_deg_params params = {alpha, (int)n_orig, R, T, estimated_coord};
  double est_edges = integrate(&exp_deg, 2, xu, xl, &params);

  return (est_edges - num_edges) * (est_edges - num_edges);
}

void NLEHelper::estimateRadialCoordinates() {
  vector<double> estimated(H.n, -1);

  FOR(i, H.n) {
    if (estimated[H.edges[i].size()] >= 0) {
      radial_coords->at(i) = estimated[H.edges[i].size()];
      continue;
    }

    // METHOD FROM PAPER
    //    radial_coords->at(i) = 2 * log(2 * n_orig * alpha * T
    //                       / (sin(M_PI * T) * H.edges[i].size() * (alpha -
    //                       0.5)));
    //    if (radial_coords->at(i) > R)
    //      radial_coords->at(i) = R;

    // TODO fix
    // NEW METHOD
    //    LOG(INFO) << "deg(" << i << ") = " << H.edges[i].size();
    radial_coords->at(i) = std::max(R + log((i + 1) / n_orig) / alpha, 0.0);
    //    LOG(INFO) << "mapped to " << radial_coords->at(i);
  }
}

#include "powerlawCommon.h"

#include <gsl/gsl_statistics_double.h>

#include "random.h"

/**
 * @author: W.M. Otte (wim@invivonmr.uu.nl); Image Sciences Institute, UMC
 * Utrecht, NL.
 * @date: 19-11-2009
 *
 * Function implementations of powerlaw scaling parameter estimation.
 *
 * ***************************************************************************
 * Method: "Power-law distributions in empirical data", Clauset et al, 2009
 * http://www.santafe.edu/~aaronc/powerlaws/
 * ***************************************************************************
 */
namespace plfit {
/**
 *
 */

void Powerlaw::SingleFit(const VectorType& V, VectorType& results, bool nosmall,
                         bool finite, float startValue, float incrementValue,
                         float endValue) {
  // determine mle-type (discrete or continuous) ...
  bool discrete = IsDiscrete(V);

  if (discrete) {
    // std::cout << "*** INFO ***: Discrete maximum likelihood estimation." <<
    // std::endl;
    MleInt(V, nosmall, finite, startValue, incrementValue, endValue, results);

  } else {
    // std::cout << "*** INFO ***: Continuous maximum likelihood estimation." <<
    // std::endl;
    MleReal(V, nosmall, finite, results);
  }
}

/**
 *
 */

void Powerlaw::BootstrapFit(const VectorType& V, VectorType& results,
                            bool nosmall, bool finiteSize, ValueType startValue,
                            ValueType incrementValue, ValueType endValue,
                            unsigned int bootstrapIterations, bool verbose) {
  // determine mle-type (discrete or continuous) ...
  bool discrete = IsDiscrete(V);

  Bootstrap(V, nosmall, finiteSize, startValue, incrementValue, endValue,
            discrete, bootstrapIterations, results, verbose);
}

/**
 * Return true of all vector elements are integers.
 */

bool Powerlaw::IsDiscrete(const VectorType& V) {
  typename VectorType::const_iterator new_end =
      std::find_if(V.begin(), V.end(), floating_point<ValueType>());
  return (V.end() == new_end);
}

/**
 * Remove duplicated values from vector and return unique values vector.
 */

void Powerlaw::Unique(const VectorType& V, VectorType& results) {
  VectorType tmp;

  std::copy(V.begin(), V.end(), std::back_inserter(tmp));

  // first sort...
  std::sort(tmp.begin(), tmp.end());

  // determine end point of unique vector elements...
  typename VectorType::iterator new_end = std::unique(tmp.begin(), tmp.end());

  // iterator...
  typename VectorType::iterator it = tmp.begin();

  // copy all unique elements in return vector...
  while (it != new_end) {
    results.push_back(*it);
    it++;
  }
}

/**
 * Remove last element of given vector.
 */

void Powerlaw::RemoveLastElement(VectorType& V) {
  V.erase(V.end() - 1, V.end());
}

/**
 * Return sorted inport vector.
 */

void Powerlaw::Sort(const VectorType& V, VectorType& W) {
  std::copy(V.begin(), V.end(), std::back_inserter(W));
  std::sort(W.begin(), W.end());
}

/**
 *
 */

void Powerlaw::KeepHigherOrEqual(VectorType& V, ValueType x) {
  typename VectorType::iterator new_end = std::remove_if(
      V.begin(), V.end(),
      std::bind(std::less<ValueType>(), std::placeholders::_1, x));

  V.erase(new_end, V.end());
}

/**
 *
 */

void Powerlaw::KeepLowerOrEqual(VectorType& V, ValueType x) {
  typename VectorType::iterator new_end = std::remove_if(
      V.begin(), V.end(),
      std::bind(std::greater_equal<ValueType>(), std::placeholders::_1, x));

  V.erase(new_end, V.end());
}

/**
 * Return incremental vector.
 */

void Powerlaw::GetIncrementVector(ValueType start, ValueType increment,
                                  ValueType end, VectorType& V) {
  int n = (end - start) / increment;

  if (n <= 0) {
    std::cerr << "*** WARNING ***: Increment vector is set to size: 0!"
              << std::endl;
  }

  for (; start <= end; start += increment) {
    V.push_back(start);
  }
}

/**
 * Return standard deviation (sqrt variance).
 */

ValueType Powerlaw::GetSD(const std::vector<ValueType>& V) {
  return gsl_stats_sd(V.data(), 1, V.size());
}

/**
 * Return cumulative sum of input vector.
 */

void Powerlaw::CumulativeSum(const VectorType& V, VectorType& W) {
  std::copy(V.begin(), V.end(), std::back_inserter(W));

  for (unsigned int i = 1; i < V.size(); i++) {
    W[i] += W[i - 1];
  }
}

/**
 * Return uniform random values vector from inputs, with similar size.
 */

void Powerlaw::GetRandomValue(const VectorType& inputs, VectorType& results) {
  for (unsigned int i = 0; i < inputs.size(); i++) {
    results.push_back(inputs[randint((int)inputs.size())]);
  }
}

/**
 * Discrete Maximum likelihood estimation.
 */

void Powerlaw::MleInt(const VectorType& x, bool nosmall, bool finiteSize,
                      ValueType startValue, ValueType increment,
                      ValueType endValue, VectorType& results) {
  VectorType vec;
  GetIncrementVector(startValue, increment, endValue, vec);
  bool finiteSizeMessagePrinted = false;

  VectorType zvec(vec.size());
  std::transform(vec.begin(), vec.end(), zvec.begin(), zeta());

  VectorType xmins;
  Unique(x, xmins);

  RemoveLastElement(xmins);

  // first and second column of data matrix...
  VectorType dat1(xmins.size());
  VectorType dat2(xmins.size());

  VectorType sorted_x;

  Sort(x, sorted_x);

  ValueType Y = 0;
  ValueType xmax = *(std::max_element(sorted_x.begin(), sorted_x.end()));

  for (unsigned int xm = 0; xm < xmins.size(); xm++) {
    ValueType xmin = xmins[xm];

    VectorType z(sorted_x);

    KeepHigherOrEqual(z, xmin);

    ValueType n = (ValueType)z.size();

    // fill L with -Inf
    VectorType L(vec.size(), -std::numeric_limits<ValueType>::infinity());

    // use copy of z, because z is used again later...
    VectorType tmp(z);
    std::transform(z.begin(), z.end(), tmp.begin(), log<ValueType>());
    ValueType slogz =
        std::accumulate(tmp.begin(), tmp.end(), static_cast<ValueType>(0));

    // xminvec = (1:xmin-1) ...
    VectorType xminvec_root(xmin - 1);
    for (int i = 0; i < xminvec_root.size(); ++i) xminvec_root[i] = 1 + i;

    for (unsigned int k = 0; k < vec.size(); k++) {
      VectorType xminvec(xminvec_root);

      ValueType exp = vec[k];

      // xminvec.^-vec(k)...
      std::transform(
          xminvec.begin(), xminvec.end(), xminvec.begin(),
          std::bind(power<ValueType>(), std::placeholders::_1, -exp));

      // sum( xminvec.^-vec(k))...
      ValueType sum = std::accumulate(xminvec.begin(), xminvec.end(),
                                      static_cast<ValueType>(0));

      // log-likelihood
      L[k] = -exp * slogz - n * std::log(zvec[k] - sum);
    }

    typename VectorType::iterator max_it = std::max_element(L.begin(), L.end());

    Y = *(max_it);
    unsigned int I = (unsigned int)std::distance(L.begin(), max_it);

    //  compute KS statistic

    ValueType exp = vec[I];

    // xmin:xmax ...
    VectorType xmin_xmax(xmax - (xmin - 1));
    for (int i = 0; i < xmin_xmax.size(); ++i) xmin_xmax[i] = xmin + i;

    // first_part = ( ( ( xmin:xmax ).^-exp ) ) ...
    VectorType first_part(xmin_xmax.size());

    std::transform(xmin_xmax.begin(), xmin_xmax.end(), first_part.begin(),
                   std::bind(power<ValueType>(), std::placeholders::_1, -exp));

    // second_part = (zvec(I) - sum( ( 1:xmin-1 ).^-exp ) ) ...

    // 1:xmin - 1 ...
    VectorType pp(xmin);
    for (int i = 0; i < pp.size(); ++i) pp[i] = 1 + i;
    RemoveLastElement(pp);

    // ( 1:xmin - 1 ) .^ - exp ...
    std::transform(pp.begin(), pp.end(), pp.begin(),
                   std::bind(power<ValueType>(), std::placeholders::_1, -exp));

    // sum( ( 1:xmin -1 ) .^ - exp ) ...
    ValueType sum_sp =
        std::accumulate(pp.begin(), pp.end(), static_cast<ValueType>(0));

    // zvec[I] - sum( ( 1:xmin -1 ) .^ - exp ) ...
    ValueType second_part = zvec[I] - sum_sp;

    // first_part / second_part ...
    std::transform(first_part.begin(), first_part.end(), first_part.begin(),
                   std::bind(std::divides<ValueType>(), std::placeholders::_1,
                             second_part));

    // fit = cumsum( first_part /. second_part ); ==
    // fit = cumsum((((xmin:xmax).^-vec(I)))./ (zvec(I) -
    // sum((1:xmin-1).^-vec(I)))) ...
    VectorType fit;
    CumulativeSum(first_part, fit);

    // normalized histogram ...
    Histogram hist(z, xmin, xmax, (xmax - xmin), true);
    VectorType histogram = hist.getHistogram();

    //  cdi = cumsum(hist(z, xmin:xmax)./n) ...
    VectorType cdi;
    CumulativeSum(histogram, cdi);

    std::transform(fit.begin(), fit.end(), cdi.begin(), fit.begin(),
                   std::minus<ValueType>());
    std::transform(fit.begin(), fit.end(), fit.begin(), abs<ValueType>());

    // dat(xm,:) = [max(abs( fit - cdi )) vec(I)] ...
    dat1[xm] = *(std::max_element(fit.begin(), fit.end()));
    dat2[xm] = vec[I];
  }

  // select the index for the minimum value of D
  // [ D, I ] = min( dat( :, 1 ) ) ...
  typename VectorType::iterator min_it =
      std::min_element(dat1.begin(), dat1.end());

  unsigned int I = (unsigned int)std::distance(dat1.begin(), min_it);
  ValueType xmin = xmins[I];

  // z = x(x>=xmin);
  // n = length(z);
  VectorType z(x);
  KeepHigherOrEqual(z, xmin);
  ValueType n = (ValueType)z.size();

  // alpha = dat( I, 2 ) ...
  ValueType alpha = dat2[I];

  // finite-size correction
  if (finiteSize) {
    alpha =
        alpha * (static_cast<ValueType>(n) - 1) / static_cast<ValueType>(n) +
        1 / static_cast<ValueType>(n);
  }

  if (!finiteSize && (n < 50) && (!finiteSizeMessagePrinted)) {
    std::cout << "*** WARNING ***: finite-size bias may be present!"
              << std::endl;
    finiteSizeMessagePrinted = true;
  }

  // L = -alpha * sum( log( z ) ) - n * log( zvec( find( vec <= alpha, 1 ,
  // 'last' ) ) - sum((1:xmin-1).^-alpha));

  // 1:xmin - 1 ...
  VectorType pp(xmin);
  for (int i = 0; i < pp.size(); ++i) pp[i] = 1 + i;
  RemoveLastElement(pp);
  // ( 1:xmin - 1 ) .^ - alpha ...
  std::transform(pp.begin(), pp.end(), pp.begin(),
                 std::bind(power<ValueType>(), std::placeholders::_1, -alpha));
  // sum( ( 1:xmin -1 ) .^ - alpha ) ...
  ValueType sum_third_part =
      std::accumulate(pp.begin(), pp.end(), static_cast<ValueType>(0));

  std::reverse(vec.begin(), vec.end());
  typename VectorType::iterator it_last = std::find_if(
      vec.begin(), vec.end(),
      std::bind(std::less_equal<ValueType>(), std::placeholders::_1, alpha));
  unsigned int index = (unsigned int)std::distance(vec.begin(), it_last) + 1;
  // n * log( zvec( find( vec <= alpha, 1 , 'last' ) ) -
  // sum((1:xmin-1).^-alpha)) ...
  ValueType sum_second_part =
      n * std::log(zvec[vec.size() - index] - sum_third_part);

  // -alpha * sum( log( z ) ) ...
  std::transform(z.begin(), z.end(), z.begin(), log<ValueType>());
  ValueType sum_first_part =
      -alpha * std::accumulate(z.begin(), z.end(), static_cast<ValueType>(0));

  results.push_back(alpha);
  results.push_back(xmin);
  results.push_back(sum_first_part - sum_second_part);
}

/**
 * Continues Maximum likelihood estimation.
 */

void Powerlaw::MleReal(const VectorType& x, bool nosmall, bool finiteSize,
                       VectorType& results) {
  VectorType xmins;
  Unique(x, xmins);
  bool finiteSizeMessagePrinted = false;

  RemoveLastElement(xmins);

  VectorType dat(xmins.size(), 0);

  VectorType sorted_x;

  Sort(x, sorted_x);

  for (unsigned int xm = 0; xm < xmins.size(); xm++) {
    VectorType z(sorted_x);

    ValueType xmin = xmins[xm];

    KeepHigherOrEqual(z, xmin);
    VectorType tmp(z);  // backup for later in fuction...

    ValueType n = (ValueType)z.size();

    // estimate alpha using direct MLE
    std::transform(
        z.begin(), z.end(), z.begin(),
        std::bind(log_div<ValueType>(), std::placeholders::_1, xmin));

    ValueType a =
        static_cast<ValueType>(n) /
        std::accumulate(z.begin(), z.end(), static_cast<ValueType>(0));

    if (nosmall) {
      if ((static_cast<ValueType>(a) - 1) / sqrt(static_cast<ValueType>(n)) >
          0.1) {
        dat.erase(dat.begin() + xm, dat.end());
        xm = (unsigned int)xmins.size() + 1;
        break;
      }
    }

    // compute KS statistic
    VectorType cx(n);
    for (int i = 0; i < cx.size(); ++i) cx[i] = i;
    std::transform(
        cx.begin(), cx.end(), cx.begin(),
        std::bind(std::divides<ValueType>(), std::placeholders::_1, n));

    // cf = xmin / z ...
    VectorType cf(tmp.size(), xmin);
    std::transform(cf.begin(), cf.end(), tmp.begin(), cf.begin(),
                   std::divides<ValueType>());

    // cf = 1 - ( xmin / z ) ^ a ...
    std::transform(
        cf.begin(), cf.end(), cf.begin(),
        std::bind(power_minus_one<ValueType>(), std::placeholders::_1, a));

    // max( abs( cf - cx ) ) ...
    std::transform(cf.begin(), cf.end(), cx.begin(), cf.begin(),
                   std::minus<ValueType>());
    dat[xm] = *(std::max_element(cf.begin(), cf.end()));
  }

  // D = min(dat) ...
  ValueType D = *(std::min_element(dat.begin(), dat.end()));

  // xmin  = xmins( find( dat <= D, 1, 'first' ) ) ...
  typename VectorType::iterator new_end = std::find_if(
      dat.begin(), dat.end(),
      std::bind(std::less_equal<ValueType>(), std::placeholders::_1, D));

  ValueType xmin = xmins[std::distance(dat.begin(), new_end)];

  // z = x( x >= xmin ) ...
  VectorType z(x);
  KeepHigherOrEqual(z, xmin);

  unsigned int n = (unsigned int)z.size();

  // alpha = 1 + n ./ sum( log(z./xmin) ) ...
  VectorType tmp(z.size(), xmin);
  std::transform(z.begin(), z.end(), tmp.begin(), z.begin(),
                 log_div<ValueType>());
  ValueType sum =
      std::accumulate(z.begin(), z.end(), static_cast<ValueType>(0));
  ValueType alpha = 1.0 + z.size() / sum;

  // finite-size correction
  if (finiteSize) {
    alpha =
        alpha * (static_cast<ValueType>(n) - 1) / static_cast<ValueType>(n) +
        1 / static_cast<ValueType>(n);
  }

  if (!finiteSize && (n < 50) && (!finiteSizeMessagePrinted)) {
    std::cout << "*** WARNING ***: finite-size bias may be present!"
              << std::endl;
    finiteSizeMessagePrinted = true;
  }

  // log-likelihood: L = n*log((alpha-1)/xmin) - alpha.*sum(log(z./xmin));
  ValueType L = n * std::log((alpha - 1) / xmin) - alpha * sum;

  results.push_back(alpha);
  results.push_back(xmin);
  results.push_back(L);
}

/**
 * Run bootstrapping of mle n-times and return Vector with indices:
 *
 * 0: average alpha
 * 1: average xmin
 * 2: average L
 *
 * 3: sd alpha
 * 4: sd xmin
 * 5: sd L
 *
 * Empty VectorType is returned when data is corrupt.
 */

void Powerlaw::Bootstrap(const VectorType& inputs, bool nosmall,
                         bool finiteSize, ValueType startValue,
                         ValueType increment, ValueType endValue, bool discrete,
                         unsigned int n, VectorType& results, bool verbose) {
  VectorType all_alpha;
  VectorType all_xmin;
  VectorType all_L;

  if (verbose)
    std::cout << "*** INFO ***: boostrapping done (%): " << std::endl;

  unsigned int successfulBootstraps = 0;

  for (unsigned int i = 0; i < n; i++) {
    if (verbose) {
      std::cout << (static_cast<ValueType>(i) / static_cast<ValueType>(n)) *
                       100;
      std::cout.flush();
      std::cout << '\r';
    }

    VectorType random_inputs;
    GetRandomValue(inputs, random_inputs);

    VectorType run;

    Mle(random_inputs, nosmall, finiteSize, startValue, increment, endValue,
        discrete, run);

    if (!run.empty()) {
      all_alpha.push_back(run[0]);
      all_xmin.push_back(run[1]);
      all_L.push_back(run[2]);
      successfulBootstraps++;
    }
  }

  if (successfulBootstraps != n) {
    ValueType p = (static_cast<ValueType>(successfulBootstraps) /
                   static_cast<ValueType>(n)) *
                  100;
    std::cerr << "*** WARNING ***: bootstrapping only ran partially "
                 "-> ("
              << p << " %)." << std::endl;
  }

  if (!all_alpha.empty()) {
    ValueType average_alpha =
        std::accumulate(all_alpha.begin(), all_alpha.end(),
                        static_cast<ValueType>(0)) /
        all_alpha.size();
    ValueType average_xmin = std::accumulate(all_xmin.begin(), all_xmin.end(),
                                             static_cast<ValueType>(0)) /
                             all_xmin.size();
    ;
    ValueType average_L =
        std::accumulate(all_L.begin(), all_L.end(), static_cast<ValueType>(0)) /
        all_L.size();
    ;

    ValueType sd_alpha = GetSD(all_alpha);
    ValueType sd_xmin = GetSD(all_xmin);
    ValueType sd_L = GetSD(all_L);

    results.push_back(average_alpha);
    results.push_back(average_xmin);
    results.push_back(average_L);

    results.push_back(sd_alpha);
    results.push_back(sd_xmin);
    results.push_back(sd_L);
  }
}

/**
 * Maximum likelihood estimation.
 *
 * If input is not correct an empty VectorType will be returned!
 */

void Powerlaw::Mle(const VectorType& inputs, bool nosmall, bool finiteSize,
                   ValueType startValue, ValueType increment,
                   ValueType endValue, bool discrete, VectorType& results) {
  // check if at all inputs are not identical ...
  if (std::max_element(inputs.begin(), inputs.end()) ==
      std::min_element(inputs.begin(), inputs.end())) {
    return;
  }

  if (discrete) {
    if (startValue > 1.0)
      MleInt(inputs, nosmall, finiteSize, startValue, increment, endValue,
             results);

    else
      std::cerr << "*** ERROR ***: start-value should be higher than 1.0!"
                << std::endl;

  } else
    MleReal(inputs, nosmall, finiteSize, results);
}
}  // namespace plfit

#pragma once

#include <gsl/gsl_sf_zeta.h>

#include <algorithm>
#include <cmath>
#include <ctime>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <numeric>
#include <vector>

//#include <boost/detail/algorithm.hpp>
//
//#include <boost/math/special_functions/zeta.hpp>
//#include <boost/math/distributions/chi_squared.hpp>
//
//#include <boost/accumulators/numeric/functional/vector.hpp>
//#include <boost/accumulators/numeric/functional/complex.hpp>
//#include <boost/accumulators/numeric/functional/valarray.hpp>
//#include <boost/accumulators/statistics/stats.hpp>
//#include <boost/accumulators/statistics/variance.hpp>
//#include <boost/accumulators/accumulators.hpp>
//#include <boost/accumulators/statistics.hpp>
//
//#include <boost/random/uniform_int.hpp>
//#include <boost/random/uniform_real.hpp>
//#include <boost/random/variate_generator.hpp>
//#include <boost/random.hpp>

/**
 * @author: W.M. Otte (wim@invivonmr.uu.nl); Image Sciences Institute, UMC
 * Utrecht, NL.
 * @date: 19-11-2009
 *
 * Function definitions of powerlaw scaling parameter estimation.
 *
 * ***************************************************************************
 * Method: "Power-law distributions in empirical data", Clauset et al, 2009
 * http://www.santafe.edu/~aaronc/powerlaws/
 * ***************************************************************************
 */
namespace plfit {

typedef double ValueType;
typedef long IntegerType;
typedef std::vector<ValueType> VectorType;

/**
 * STL Predicate: find floating number.
 */
template <class T>
struct floating_point : public std::unary_function<T, bool> {
  bool operator()(const T& x) const { return !(floor(x) == x); }
};

/**
 * STL helper: print number.
 */
template <class T>
struct print : public std::unary_function<T, void> {
  void operator()(const T& x) const { std::cout << x << std::endl; }
};

/**
 * STL helper: log( x / y ).
 */
template <class T>
struct log_div : public std::binary_function<T, T, T> {
  T operator()(const T& x, const T& y) const { return std::log(x / y); }
};

/**
 * STL helper: log( x ).
 */
template <class T>
struct log : public std::unary_function<T, T> {
  T operator()(const T& x) const { return std::log(x); }
};

/**
 * STL helper: 1 - x^N.
 */
template <class T>
struct power_minus_one : public std::binary_function<T, T, T> {
  T operator()(const T& x, const T& N) const { return 1.0f - pow(x, N); }
};

/**
 * STL helper: x^N.
 */
template <class T>
struct power : public std::binary_function<T, T, T> {
  T operator()(const T& x, const T& N) const { return pow(x, N); }
};

/**
 * Boost library zeta function.
 */
//	template< class T >
//	struct zeta: public std::unary_function< T, T >
//	{
//		T operator()( const T& x ) const
//		{
//      return boost::math::zeta< T >( x );
//		}
//	};
//   template<>
struct zeta : public std::unary_function<double, double> {
  double operator()(const double& x) const { return gsl_sf_zeta(x); }
};

/**
 * abs.
 */
template <class T>
struct abs : public std::unary_function<T, T> {
  T operator()(const T& x) const { return static_cast<T>(std::fabs(x)); }
};

// *********************************************************
// POWERLAW
// *********************************************************

/**
 *
 */
class Powerlaw {
 public:
  //		typedef gsl_rng random_number_type;
  //		typedef boost::uniform_real< ValueType > real_distribution_type;
  //		typedef boost::uniform_int< IntegerType > int_distribution_type;
  //		typedef boost::variate_generator< random_number_type&,
  // real_distribution_type > real_generator_type; 		typedef
  // boost::variate_generator< random_number_type&, int_distribution_type >
  // int_generator_type;

  /**
   *
   */
  static void SingleFit(const VectorType& input, VectorType& results,
                        bool nosmall, bool finiteSize, float startValue,
                        float incrementValue, float endValue);

  /**
   *
   */
  static void BootstrapFit(const VectorType& input, VectorType& results,
                           bool nosmall, bool finiteSize, ValueType startValue,
                           ValueType incrementValue, ValueType endValue,
                           unsigned int bootstrapIterations, bool verbose);

 protected:
  /**
   *
   */
  static bool IsDiscrete(const VectorType& V);

  /**
   *
   */
  static void Unique(const VectorType& V, VectorType& results);

  /**
   *
   */
  static void RemoveLastElement(VectorType& V);

  /**
   *
   */
  static void Sort(const VectorType& V, VectorType& W);

  /**
   *
   */
  static void KeepLowerOrEqual(VectorType& V, ValueType x);

  /**
   *
   */
  static void KeepHigherOrEqual(VectorType& V, ValueType x);

  /**
   *
   */
  static void GetIncrementVector(ValueType s, ValueType i, ValueType e,
                                 VectorType& V);

  /**
   *
   */
  static ValueType GetSD(const std::vector<ValueType>& V);

  /**
   *
   */
  static void CumulativeSum(const VectorType& V, VectorType& W);

  /**
   *
   */
  static void GetRandomValue(const VectorType& inputs, VectorType& results);

  /**
   *
   */
  static void MleInt(const VectorType& x, bool nosmall, bool finiteSize,
                     ValueType startValue, ValueType increment,
                     ValueType endValue, VectorType& results);

  /**
   *
   */
  static void MleReal(const VectorType& x, bool nosmall, bool finiteSize,
                      VectorType& results);

  /**
   *
   */
  static void Bootstrap(const VectorType& inputs, bool nosmall, bool finiteSize,
                        ValueType startValue, ValueType increment,
                        ValueType endValue, bool discrete, unsigned int n,
                        VectorType& results, bool verbose);

  /**
   *
   */
  static void Mle(const VectorType& inputs, bool nosmall, bool finiteSize,
                  ValueType startValue, ValueType increment, ValueType endValue,
                  bool discrete, VectorType& results);
};

/**
 * Simple histogram implementation.
 */

class Histogram {
 private:
  std::vector<ValueType> histogram;

 public:
  /**
   * Construct histogram.
   */
  Histogram(const std::vector<ValueType>& data, ValueType low, ValueType high,
            unsigned int bins, bool normalize) {
    histogram.assign(bins + 1, 0);  // plus outlier bin...
    ValueType width = (high - low) / static_cast<ValueType>(bins);

    for (unsigned int i = 0; i < data.size(); i++) {
      unsigned int bin = static_cast<unsigned int>((data[i] - low) / width);

      if (bin < bins)
        histogram[bin]++;

      else if (bin >= bins)  // insert outliers in highest bin...
        histogram[bins]++;
    }

    if (normalize)
      std::transform(histogram.begin(), histogram.end(), histogram.begin(),
                     std::bind(std::divides<ValueType>(), std::placeholders::_1,
                               data.size()));
  }

  /**
   * Return histogram.
   */
  std::vector<ValueType> getHistogram() { return histogram; }
};

}  // namespace plfit

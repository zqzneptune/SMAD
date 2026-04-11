#ifndef SPC_STATS_FNS_HPP_
#define SPC_STATS_FNS_HPP_

#include <valarray>
#include <vector>
#include <cmath>
#include <boost/foreach.hpp>

namespace saint_spc {

typedef unsigned short Count_t;

const unsigned log_factorial_max = 10000;
extern const double log_factorial_table[log_factorial_max + 1];

inline double log_factorial(const unsigned short n){
	return (n > log_factorial_max)? lgamma(n+1) : log_factorial_table[n];
}

inline double GP_log_pmf(const unsigned short x, const double lambda1, const double lambda2) {
	if(x == 0) return -lambda1;
	double tmp = lambda1 + x * lambda2;
	if (tmp <= 0) return -100; // safety
	return log(lambda1) + (x - 1) * log(tmp) - tmp - log_factorial(x);
}

// generalized poisson parameterized by mean and lambda2
inline double GP_log_pmf1(const Count_t k, const double mean, const double lambda2) {
	const double lambda1 = mean * (1 - lambda2);
	return GP_log_pmf(k, lambda1, lambda2);
}

inline double pois_log_pmf(Count_t k, const double lambda) {
	if(k == 0) return -lambda;
	return k * log(lambda) - lambda - log_factorial(k);
}

template<typename T>
double mean(const std::valarray<T>& v) {
	if (v.size() == 0) return 0.0;
	return static_cast<double>(v.sum())/v.size();
}

template<typename T>
double var1(const std::valarray<T>& v) {
	if(v.size()<=1)
		return 0.0;
	std::valarray<double> tmp(v.size());
	for(size_t i=0; i<v.size();i++)
		tmp[i] = v[i];
	tmp -= mean(v);
	return (tmp*tmp).sum()/(v.size()-1.0);
}

template<typename T>
double var1(const std::vector<T>& v) {
	if (v.size() == 0) return 0.0;
	std::valarray<double> tmp(v.size());
	for(size_t i=0; i<v.size();i++)
		tmp[i] = v[i];
	return var1(tmp);
}

} // namespace saint_spc

#endif // SPC_STATS_FNS_HPP_

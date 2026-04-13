#ifndef INT_STATS_FNS_HPP_
#define INT_STATS_FNS_HPP_

// [[Rcpp::depends(BH)]]
#include <boost/math/distributions/normal.hpp>

namespace saint_int {

// log(√2π)
// log CDF of normal by approximation
inline double normal_cdf(const double x, const double mu=0, const double sigma=1) {
	using namespace boost::math;
	return cdf(normal(mu, sigma), x);
}

const double log_sqrt_2_pi = .918938533204672741780329736405617639861397;
inline double normal_log_pdf(const double x, const double mu=0, const double sigma=1) {
	const double z = (x-mu)/sigma;
	return -.5*z*z - log(sigma) - log_sqrt_2_pi;
}
inline double normal_log_cdf(const double x, const double mu=0, const double sigma=1) {
	using namespace boost::math;
	return log(cdf(normal(mu, sigma), x));
}

template<typename T>
double mean(const std::valarray<T>& v) {
	return static_cast<double>(v.sum())/v.size();
}
template<typename T>
double mean(const std::vector<T>& v) {
	std::valarray<double> tmp(v.size());
	for(size_t i=0; i<v.size();i++)
		tmp[i] = v[i];
	return mean(tmp);
}

// the function returns sample variance if it can be calculated, and return zero otherwise
template<typename T>
double var1(const std::valarray<T>& v) {
	if(v.size()<=1)
		throw std::runtime_error("sample variance from n<2");
	std::valarray<double> tmp(v.size());
	for(size_t i=0; i<v.size();i++)
		tmp[i] = v[i];
	tmp -= mean(v);
	return (tmp*tmp).sum()/(v.size()-1.0);
}
template<typename T>
double var_MLE(const std::valarray<T>& v) {
	if(v.size()<=1)
		return 0.0;
	std::valarray<double> tmp(v.size());
	for(size_t i=0; i<v.size();i++)
		tmp[i] = v[i];
	tmp -= mean(v);
	return (tmp*tmp).sum()/v.size();
}

template<typename T>
double var1(const std::vector<T>& v) {
	std::valarray<double> tmp(v.size());
	for(size_t i=0; i<v.size();i++)
		tmp[i] = v[i];
	return var1(tmp);
}

} // namespace saint_int

#endif

#ifndef SPC_NLOPT_WRAP_HPP_
#define SPC_NLOPT_WRAP_HPP_

#include <nlopt.hpp>
#include <functional>

namespace saint_spc {

static std::vector<double> dummy_vector;

// this function does void pointer casting
inline double obj_f_wrapper(const std::vector<double>& x, std::vector<double>& grad, void* f_ptr){
	const auto& bind_f = *static_cast<std::function<double(const std::vector<double>&, std::vector<double>&)>* >(f_ptr);
	return bind_f(x,grad);
}

template<typename Func_t, typename this_type, typename... Params>
double nlopt_wrap(nlopt::opt& opt, std::vector<double>& x, Func_t func, this_type* this_ptr, const Params&... params) {
	using namespace std::placeholders;
	std::function<double(const std::vector<double>&, std::vector<double>&)> bind_f = std::bind(func, std::cref(*this_ptr),_1, _2, std::cref(params)...);
	opt.set_max_objective(obj_f_wrapper, &bind_f);
	double maxf;
	opt.optimize(x, maxf);
	return maxf;
}

} // namespace saint_spc

#endif // SPC_NLOPT_WRAP_HPP_

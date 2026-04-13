/*
 * spc_stats.cpp
 * Corrected CPC statistical model to match original SAINTexpress logic.
 */

#include <iostream>
#include <iomanip>
#include <string>
#include <valarray>
#include <vector>
// [[Rcpp::depends(BH)]]
#include <boost/array.hpp>
#include <algorithm>
#include <numeric>
#include <boost/algorithm/string.hpp>
#include <nlopt.hpp>
#include "SpcGlobals.hpp"
#include "SpcStats_fns.hpp"
#include "SpcNloptWrap.hpp"

namespace saint_spc {

using namespace std;

namespace saint {
	bool iszero(Count_t a) {return a==0;}
}

Model_data statModel(const Options& opts, const vector<vector<size_t> >& p2p_mapping, const vector<string>& ubait, const Countmat& test_mat_DATA, const Countmat& ctrl_mat_DATA, const vector<size_t>& ip_idx_to_bait_no, const size_t nprey, const size_t nbait) {

	// bait index to ip index
	vector<vector<size_t> > bait_no_to_ip_idxes(ubait.size());
	for(size_t i = 0; i < ubait.size(); ++i)
		for(size_t j = 0; j < ip_idx_to_bait_no.size(); ++j)
			if( ip_idx_to_bait_no[j] == i )
				bait_no_to_ip_idxes[i].push_back(j);

	vector<unsigned char> n_rep_vec(nbait);
	transform(bait_no_to_ip_idxes.begin(), bait_no_to_ip_idxes.end(),
			  n_rep_vec.begin(),
			  container_size<vector<size_t> >);

	const size_t n_ctrl_ip = ctrl_mat_DATA.n_cols;
	vector<double> eta(nprey), d(nprey);
	
	int L_user_input = opts.L;
	if(L_user_input <= 0) {
		throw runtime_error("L must be a positive integer.");
	}
	int L = min<size_t>(n_ctrl_ip, L_user_input);

	for(size_t idx = 0; idx < ctrl_mat_DATA.n_rows; idx++){
		valarray<Count_t> ctrl_va = ctrl_mat_DATA(idx, __);
		vector<Count_t> ctrl(begin(ctrl_va), end(ctrl_va));
		sort(ctrl.begin(), ctrl.end(), greater<Count_t>());
		
		Count_t a = accumulate(ctrl.begin(), ctrl.begin() + min<size_t>(ctrl.size(), L), 0);
		double tmp = static_cast<double>(a)/L;
		if(tmp < 0.1) tmp = 0.1;
		eta[idx] = tmp;
	}

	for(size_t i=0; i < nprey; ++i) {
		mu_init: 
		double mu = 5.0 * eta[i];
		d[i] = mu - eta[i];
	}

	// matrix of valarrays on original data
	Fastmat<valarray<Count_t> > test_mat(nprey, ubait.size());
	for(size_t i = 0; i < test_mat.n_rows; ++i)
		for(size_t j = 0; j < test_mat.n_cols; ++j){
			test_mat(i,j).resize(bait_no_to_ip_idxes.at(j).size());
			for(size_t k = 0; k < bait_no_to_ip_idxes.at(j).size(); ++k) {
				const auto ip_idx = bait_no_to_ip_idxes.at(j).at(k);
				test_mat(i,j)[k] = test_mat_DATA(i,ip_idx);
			}
		}

	Countmat test_mat_DATA1(test_mat_DATA);
	Fastmat<valarray<Count_t> > test_mat1(nprey, ubait.size());
	for(size_t i = 0; i < test_mat1.n_rows; ++i)
		for(size_t j = 0; j < test_mat1.n_cols; ++j){
			test_mat1(i,j).resize(bait_no_to_ip_idxes.at(j).size());
			for(size_t k = 0; k < bait_no_to_ip_idxes.at(j).size(); ++k) {
				const auto ip_idx = bait_no_to_ip_idxes.at(j).at(k);
				test_mat1(i,j)[k] = test_mat_DATA1(i,ip_idx);
			}
		}

	Countmat test_SS1(nprey, ubait.size());
	for(size_t i = 0; i < test_SS1.n_rows; ++i)
		for(size_t j = 0; j < test_SS1.n_cols; ++j)
			test_SS1(i, j) = test_mat1(i, j).sum();

	// initialize Z
	Fastmat<valarray1<bool> > Z(nprey, ubait.size());
	for(size_t i=0; i < nprey; ++i)
		for(size_t j=0; j < nbait; ++j)
			Z(i,j).resize(n_rep_vec.at(j));

	for(size_t i=0; i < nprey; ++i)
		for(size_t j=0; j < nbait; ++j)
			for(size_t rep = 0; rep < bait_no_to_ip_idxes.at(j).size(); ++rep)
				Z(i,j)[rep] = (test_mat1(i,j)[rep] > (2 * eta[i]));

	// calculate vector of control means
	vector<Count_t> ctrl_mean(nprey);
	for (unsigned i=0; i < nprey; ++i) {
		valarray<Count_t> row = ctrl_mat_DATA(i, __);
		ctrl_mean.at(i) = static_cast<double>(row.sum()) / row.size();
	}

	valarray<bool> apply_MRF(nprey);
	for (size_t i = 0; i < nprey; i++) {
		size_t sum = 0;
		for(size_t bait_idx = 0; bait_idx < nbait; bait_idx++)
			for(size_t rep=0; rep < Z(i,bait_idx).size(); rep++)
				sum += Z(i,bait_idx)[rep];
		if(sum < test_mat_DATA.n_cols * opts.f)
			apply_MRF[i] = true;
	}

	// lambda2 for Z=0
	valarray<double> lambda2_false(0.1, nprey);
	for(size_t i = 0; i < nprey; ++i) {
		valarray<Count_t> ctrl_va = ctrl_mat_DATA(i, __);
		vector<Count_t> ctrl(begin(ctrl_va), end(ctrl_va));
		sort(ctrl.begin(), ctrl.end(), greater<Count_t>());
		
		vector<Count_t> sub_ctrl(ctrl.begin(), ctrl.begin() + min<size_t>(ctrl.size(), L));
		double variance = var1(sub_ctrl);
		if(variance > eta[i])
			lambda2_false[i] = 1 - sqrt(eta[i]/variance);
	}

	valarray<double> lambda2_true(0.1, nprey);
	double beta0 = 0.0, beta1 = 0, gamma = 0.;

	return {
		nprey, nbait,
		test_mat_DATA,
		test_mat,
		test_mat1,
		test_SS1,
		n_rep_vec,
		p2p_mapping,
		ctrl_mean,
		static_cast<size_t>(n_ctrl_ip),
		apply_MRF,
		Z,
		beta0, beta1, gamma,
		eta, d, lambda2_true, lambda2_false
	};
}

double Model_data::llik_MRF_gamma_0( const std::vector<double>& x, std::vector<double>& /*grad*/) const {
	const double beta0_ = 0;
	double beta1_ = x[0];
	double loglik = 0;
	double MRFtrue = exp(beta1_);
	const double MRFfalse = exp(beta0_);
	for(size_t i=0; i < nprey; ++i)
		for(size_t j=0; j < nbait; ++j) {
			const auto& y = test_mat1(i,j);
			if(all_of(begin(y), end(y), saint::iszero)) continue;
			const auto& k = Z(i,j);
			const auto m = n_rep_vec[j];
			double prod = 1;
			for(unsigned rep = 0; rep < m; rep++) {
				prod *= (k[rep] ? MRFtrue : MRFfalse)/ (MRFtrue + MRFfalse);
			}
			loglik += log(prod)/m;
		}
	return loglik;
}

double Model_data::llik_MRF( const std::vector<double>& x, std::vector<double>& /*grad*/, const Fastmat<double>& gsum_mat) const {
	const double beta0_ = 0;
	double beta1_ = x[0], gamma_ = x[1];
	const double MRFfalse = exp(beta0_);
	double loglik = 0;
	for(size_t i=0; i < nprey; ++i)
		for(size_t j=0; j < nbait; ++j) {
			const auto& y = test_mat1(i,j);
			if(all_of(begin(y), end(y), saint::iszero)) continue;
			const auto& k = Z(i,j);
			const auto m = n_rep_vec[j];
			double MRFtrue = exp(beta1_ + gamma_ * gsum_mat(i,j));
			double prod = 1;
			for(unsigned rep = 0; rep < m; rep++)
				prod *= (k[rep] ? MRFtrue : MRFfalse)/ (MRFtrue + MRFfalse);
			loglik += log(prod)/m;
		}
	return loglik;
}

double Model_data::llikelihood() const {
	double loglik = 0;
	for(size_t i=0; i < nprey; ++i) {
		for(size_t j=0; j < nbait; ++j) {
			const auto& y = test_mat1(i,j);
			if(all_of(begin(y), end(y), saint::iszero)) continue;
			const auto& k = Z(i,j);
			const auto m = n_rep_vec[j];
			double gsum = 0;
			BOOST_FOREACH(const auto l , p2p_mapping[i])
				gsum += Z(l, j).mean();
			double MRFtrue = exp(beta1 + gamma * gsum);
			double MRFfalse = exp(beta0);
			double prod = 1;
			for(unsigned rep = 0; rep < m; rep++){
				prod *= (k[rep] ? MRFtrue : MRFfalse ) / (MRFtrue + MRFfalse) * 
						exp(k[rep] ? GP_log_pmf1( min<Count_t>(y[rep], static_cast<Count_t>(eta[i]+d[i])), (eta[i]+d[i]), lambda2_true[i]) : 
									 GP_log_pmf1( max<Count_t>(y[rep], static_cast<Count_t>(eta[i])), eta[i], lambda2_false[i]));
			}
			loglik += log(prod)/m;
		}
	}
	return loglik;
}

double Model_data::loglikelihood_Z(const size_t i, const size_t j, const size_t rep, const Fastmat<vector<boost::array<double, 2> > >& pre_calc_loglik) const {
	const auto& k = Z(i,j);
	double gsum = 0;
	if(gamma!=0){
		BOOST_FOREACH(const auto l , p2p_mapping[i])
			gsum += Z(l, j).mean();
	}
	double logMRFtrue = (beta1 + gamma * gsum);
	double logMRFfalse = (beta0);
	const boost::array<double, 2>& pcl = pre_calc_loglik(i,j)[rep];
	return k[rep]? logMRFtrue + pcl[0] : logMRFfalse + pcl[1];
}

void Model_data::icm_Z() {
	Fastmat<vector<boost::array<double, 2> > > pre_calc_loglik(nprey, nbait);
	for(size_t i=0; i < nprey; ++i)
		for(size_t j=0; j < nbait; ++j) {
			const auto& y = test_mat1(i,j);
			const auto m = n_rep_vec[j];
			pre_calc_loglik(i,j).resize(m);
			for(unsigned rep = 0; rep < m; rep++)
				pre_calc_loglik(i,j)[rep] = boost::array<double, 2>{{
						GP_log_pmf1(min<Count_t>(y[rep], static_cast<Count_t>(eta[i]+d[i])), (eta[i]+d[i]), lambda2_true[i]),
						GP_log_pmf1(max<Count_t>(y[rep], static_cast<Count_t>(eta[i])), eta[i], lambda2_false[i])}};
		}
	for(unsigned iter_Z=0; iter_Z < 1; ++iter_Z)
		for(size_t i=0; i < nprey; ++i)
			for(size_t j=0; j < nbait; ++j) {
				const auto m = n_rep_vec[j];
				for(unsigned rep = 0; rep < m; rep++){
					double first = loglikelihood_Z(i,j,rep,pre_calc_loglik);
					Z(i,j)[rep] = !Z(i,j)[rep];
					double second = loglikelihood_Z(i,j,rep,pre_calc_loglik);
					if(first > second) Z(i,j)[rep] = !Z(i,j)[rep];
				}
			}
}

void Model_data::wrt_MRF_gamma_0() {
	nlopt::opt opt(nlopt::LN_COBYLA, 1);
	opt.set_lower_bounds({-15});
	opt.set_upper_bounds({ 15});
	opt.set_ftol_abs(1e-4);
	opt.set_maxeval(1e2);
	std::vector<double> x{beta1};
	double oldbeta1 = beta1;
	static std::vector<double> dummy_vector;
	double oldf = llik_MRF_gamma_0(x, dummy_vector);
	double maxf = nlopt_wrap(opt, x, &Model_data::llik_MRF_gamma_0, this);
	beta1 = x[0];
	if(maxf < oldf)
		beta1 = oldbeta1;
}

void Model_data::wrt_MRF() {
	Fastmat<double> gsum_mat(nprey, nbait);
	for(size_t i=0; i < nprey; ++i)
		for(size_t j=0; j < nbait; ++j)
			BOOST_FOREACH(const auto l , p2p_mapping[i])
				gsum_mat(i,j) += Z(l, j).mean();
	nlopt::opt opt(nlopt::LN_COBYLA, 2);
	opt.set_lower_bounds({-15, 0});
	opt.set_upper_bounds({ 15, 10});
	opt.set_ftol_abs(1e-4);
	opt.set_maxeval(1e2);
	std::vector<double> x{beta1, gamma};
	double oldbeta1 = beta1, oldgamma = gamma;
	static std::vector<double> dummy_vector;
	double oldf = llik_MRF(x, dummy_vector, gsum_mat);
	double maxf = nlopt_wrap(opt, x, &Model_data::llik_MRF, this, gsum_mat);
	beta1 = x[0], gamma = x[1];
	if(maxf <= llik_MRF({x[0],0},dummy_vector, gsum_mat))
		gamma = 0;
	if(maxf < oldf)
		beta1 = oldbeta1, gamma = oldgamma;
}

void Model_data::wrt_d() {
	for(size_t i=0; i < nprey; ++i) {
		Count_t sum = 0;
		unsigned sum_count = 0;
		for(size_t j=0; j < nbait; ++j) {
			const auto& k = Z(i,j);
			const auto& y = test_mat1(i,j);
			const auto m = n_rep_vec[j];
			for(unsigned rep = 0; rep < m; rep++)
				if(k[rep]){
					sum += y[rep];
					sum_count++;
				}
		}
		if(sum_count>0)	d[i] = static_cast<double>(sum)/sum_count - eta[i];
		else d[i] = 4 * eta[i];
	}
}

boost::array<Fastmat<double>, 3> Model_data::calculateScore(const Options& opts) const {
	Fastmat<double> average_score(nprey, nbait);
	Fastmat<double> min_log_odds_score(nprey, nbait);
	Fastmat<double> maximum_score(nprey, nbait);
	for(size_t i=0; i < nprey; ++i)
		for(size_t j=0; j < nbait; ++j) {
			const auto& y = test_mat1(i,j);
			const auto m = n_rep_vec[j];
			double gsum = 0;
			BOOST_FOREACH(const auto l , p2p_mapping[i])
				gsum += Z(l, j).mean();
			double MRFtrue = exp(beta1 + gamma * gsum);
			double MRFfalse = exp(beta0);
			average_score(i, j) = 0;
			min_log_odds_score(i, j) = 0;
			vector<double> tmp_scores(m),tmp_odds_scores(m);
			for(unsigned rep = 0; rep < m; ++rep) {
				double tmp_mean = 0;
				if(y[rep] >= 5) tmp_mean = min<double>(eta[i] + d[i], 5 * eta[i]);
 				else tmp_mean = (eta[i] + d[i]);

				double unnorm_score_true = MRFtrue;
				double unnorm_score_false = MRFfalse * exp(GP_log_pmf1(max<Count_t>(y[rep], static_cast<Count_t>(eta[i])), eta[i], lambda2_false[i]) - 
														  GP_log_pmf1(min<Count_t>(y[rep], static_cast<Count_t>(tmp_mean)), tmp_mean, lambda2_true[i]));
				
				if(y[rep] <= 1 || static_cast<double>(y[rep]) <= ctrl_mean.at(i)){
					tmp_scores.at(rep) = 0;
				}else{
					tmp_scores.at(rep) = unnorm_score_true / (unnorm_score_true + unnorm_score_false);
				}
				tmp_odds_scores.at(rep)=log(unnorm_score_true)-log(unnorm_score_false);
			}
			const size_t max_rep = opts.R;
			std::sort(tmp_scores.begin(), tmp_scores.end(), greater<double>());
			std::sort(tmp_odds_scores.begin(), tmp_odds_scores.end(), greater<double>());
			if(m > max_rep){
				average_score(i, j) = std::accumulate(tmp_scores.begin(), tmp_scores.begin() + max_rep, 0.0 ) / ((double) max_rep);
			}else{
				average_score(i, j) = std::accumulate(tmp_scores.begin(), tmp_scores.end(), 0.0 ) / ((double) m);
			}
			min_log_odds_score(i, j) = tmp_odds_scores.back();
			maximum_score(i, j) = tmp_scores.at(0);

			if(y.max()==1)
				average_score(i,j) = maximum_score(i,j) = 0;
		}
	return boost::array<Fastmat<double>, 3>{{average_score, maximum_score, min_log_odds_score}};
}

} // namespace saint_spc

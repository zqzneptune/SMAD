#ifndef INT_GLOBALS_HPP_
#define INT_GLOBALS_HPP_

#include <Rcpp.h>

#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <set>
#include <cmath>

#include <boost/foreach.hpp>
#include <boost/array.hpp>
#include "../saint_common/PreyClass.hpp"
#include "../saint_common/BaitClass.hpp"
#include "../saint_common/InterClass.hpp"
#include "../saint_common/UIClass.hpp"

namespace saint_int {

typedef unsigned Count_t;
class Quant_t;

namespace saint {
	inline bool isnan(const Quant_t);
}

class Quant_t {
	friend Quant_t exp(Quant_t);
	friend bool saint::isnan(const Quant_t x);
 	double quant;
public:
 	Quant_t(double x) : quant(x){ }
 	Quant_t() : quant(0){ }
 	inline Quant_t& operator= (double a) {
		quant = a;
		return *this;
	}
 	inline operator double& () {
 		return quant;
 	}
 	inline operator const double& () const {
 		return quant;
 	}
	inline bool operator<(Quant_t r) const {
		return quant < r.quant;
	}
};

} // namespace saint_int

// Fastmat needs Count_t and Quant_t defined before inclusion
#include "IntFastmat.hpp"
#include "IntStats_fns.hpp"

namespace saint_int {

using namespace std;

const unsigned iter_Z_d = 10;

namespace saint {
	const double nan = std::numeric_limits<double>::quiet_NaN();
	inline bool isnan(const Quant_t x) { return std::isnan(x.quant); }
}

inline Quant_t exp(Quant_t a){
	if(saint::isnan(a)) return saint::nan;
	return std::exp(a.quant);
}

vector<string> splitString(const string & input_txt);
void getFileDimensions(string inputFile, int &nr, int &nc);

deque<InterClass> parseInterFile(const string & inputFile, size_t &ninter);

void mapRowCol( deque<InterClass> &IDATA, const deque<PreyClass> &PDATA, const deque<BaitClass> &BDATA, const map<string, BaitClass>&);

std::set<string> parsePreyFile( deque<PreyClass> &PDATA, const string & inputFile, size_t &nprey);

map<string, BaitClass> parseBaitFile( deque<BaitClass> &BDATA, const string & inputFile, size_t &nip);
void sortBaitData( deque<BaitClass> &BDATA );

void createList( deque<UIClass> &UIDATA, const deque<InterClass> &IDATA, const deque<BaitClass> &BDATA, size_t &nuinter );

int get_nexpr( deque<BaitClass> &BDATA );
int get_nctrl( const deque<BaitClass> &BDATA );

inline std::string detect_line_ending(const string& file_name) {
	std::ifstream infile(file_name);
	char tmp;
	for(int i = 0;i < 1e3; i++){
		infile.get(tmp);
		if(tmp == '\r') {
			infile.get(tmp);
			if(tmp == '\n' )
				return "\r\n";
			return "\r";
		}
		if(tmp == '\n')
			return "\n";
	}
	return "";
}

template<typename T>
size_t container_size(const T& container) {
	return container.size();
}

inline istream& skip_line(istream& is, const string& eol, const unsigned times = 1){
	for(unsigned i=0; i<times; i++)
		is.ignore(std::numeric_limits<std::streamsize>::max(), eol[eol.size()-1]);
	return is;
}

inline double quant_log_pdf(Quant_t x, const double mu=0, const double sigma=1) {
	if(saint::isnan(x)) {
		throw std::runtime_error("isnan in quant_log_pdf");
	} else {
		if(abs(x-mu)/sigma>5)
			x=(signbit(static_cast<double>(x))?-1:1)*5*sigma+mu;
		return normal_log_pdf(x,mu,sigma);
	}
}

struct Options {
	double f;
	unsigned R;
	int L;
};

struct Model_data {
	// data
	const size_t nprey;
	const size_t nbait;

	const Quantmat test_mat_DATA;
	const Fastmat<valarray<Quant_t> > test_mat;
	const Fastmat<valarray<Quant_t> > test_mat1;
	const vector<unsigned char> n_rep_vec;
	const vector<vector<size_t> > p2p_mapping;
	const double t;
	const size_t n_ctrl_ip;

	// parameters
	Fastmat<valarray1<bool> > Z;

	double beta0, beta1, gamma;

	const vector<double> eta;
	vector<double> d;
	const valarray<double> sd_true;
	const valarray<double> sd_false;


	void print_data_summary(const deque<PreyClass> &PDATA,const vector<string> &ubait, const deque<UIClass> &UIDATA) const;

	double llikelihood() const;
	double loglikelihood_Z(const size_t i, const size_t j, const size_t rep, const Fastmat<vector<boost::array<double, 2> > >& pre_calc_loglik) const;
	boost::array<Fastmat<double>, 3> calculateScore(const Options& opts) const;
	vector<double> get_MRF_parameters() const {
		return {beta0,beta1,gamma};
	}
	void set_MRF_parameters(double beta0_, double beta1_, double gamma_){
		beta0=beta0_, beta1=beta1_, gamma=gamma_;
	}
	void set_MRF_parameters(const vector<double>& a){
		if(a.size()!=3)
			throw runtime_error("3 MRF parameters needed");
		beta0=a[0], beta1=a[1], gamma=a[2];
	}
	void print_MRF_parameters(ostream& out = Rcpp::Rcout) const {
		out << "";
	}
	void wrt_MRF();
	void wrt_MRF_gamma_0();
	void wrt_d();
	void icm_Z();
private:
	Fastmat<vector<boost::array<double, 2> > > precalculate_densities() const;
	double llik_MRF_gamma_0( const std::vector<double> &x, std::vector<double>& /*grad*/) const;
	double llik_MRF( const std::vector<double>& x, std::vector<double>& /*grad*/, const Fastmat<double>& gsum) const;
	double llik_gamma( const std::vector<double>& x, std::vector<double>& /*grad*/, const Fastmat<double>& gsum, const Fastmat<vector<boost::array<double, 2> > >& log_densities) const;
};

} // namespace saint_int

#endif /* INT_GLOBALS_HPP_ */

#ifndef SPC_GLOBALS_HPP_
#define SPC_GLOBALS_HPP_

#include <Rcpp.h>

#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <set>

#include <array>
#include "../saint_common/PreyClass.hpp"
#include "../saint_common/BaitClass.hpp"
#include "../saint_common/InterClass.hpp"
#include "../saint_common/UIClass.hpp"

namespace saint_spc {

typedef unsigned short Count_t;

} // namespace saint_spc

// Fastmat needs Count_t defined before inclusion
#include "SpcFastmat.hpp"

namespace saint_spc {

using namespace std;

double str2dbl(string ch);

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
		if(tmp == '\r') {		// Mac or Windows
			infile.get(tmp);
			if(tmp == '\n' )	// Windows
				return "\r\n";
			return "\r";		// Mac
		}
		if(tmp == '\n')			// Unix
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

struct Options {
	double f;
	unsigned R;
	int L;
};

struct Model_data {
	// data
	const size_t nprey, nbait;
	const Countmat test_mat_DATA;
	const Fastmat<valarray<Count_t> > test_mat;
	const Fastmat<valarray<Count_t> > test_mat1;
	const Countmat test_SS;
	const vector<unsigned char> n_rep_vec;
	const vector<vector<size_t> > p2p_mapping;
	const vector<Count_t> ctrl_mean;
	const size_t n_ctrl_ip;
	const valarray<bool> apply_MRF;

	// parameters
	Fastmat<valarray1<bool> > Z;
	double beta0, beta1, gamma;

	const vector<double> eta;
	vector<double> d;
	const valarray<double> lambda2_true, lambda2_false;

	void print_data_summary(const deque<PreyClass> &PDATA,const vector<string> &ubait, const deque<UIClass> &UIDATA) const;

	double llikelihood() const;
	double loglikelihood_Z(const size_t i, const size_t j, const size_t rep, const Fastmat<vector<std::array<double, 2> > >& pre_calc_loglik) const;
	std::array<Fastmat<double>, 3> calculateScore(const Options& opts) const;
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
	Fastmat<vector<std::array<double, 2> > > precalculate_poisson_pmfs() const;
	double llik_MRF_gamma_0( const std::vector<double> &x, std::vector<double>& /*grad*/) const;
	double llik_MRF( const std::vector<double>& x, std::vector<double>& /*grad*/, const Fastmat<double>& gsum) const;
	double llik_gamma( const std::vector<double>& x, std::vector<double>& /*grad*/, const Fastmat<double>& gsum, const Fastmat<vector<std::array<double, 2> > >& poisson_pmfs) const;
};

} // namespace saint_spc

#endif /* SPC_GLOBALS_HPP_ */

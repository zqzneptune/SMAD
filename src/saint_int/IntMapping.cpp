/*
 * int_Mapping.cpp
 * Adapted from SAINTexpress SAINT-MRF-int/Mapping.cpp for Rcpp integration
 */

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <algorithm>
#include <set>
#include <boost/algorithm/string.hpp>
#include "IntGlobals.hpp"
#include "IntStats_fns.hpp"

namespace saint_int {

using namespace std;

vector<string> uniqueBait(const deque<BaitClass> &BDATA, size_t &nbait ) {
	set<string> ubait_set;
	vector<string> ubait_vec;
	BOOST_FOREACH(const BaitClass& bdata , BDATA)
		if( !bdata.get_isCtrl() ){
			ubait_vec.push_back(bdata.get_baitId());
			ubait_set.insert(bdata.get_baitId());
		}
	{
		auto tmp_iter = std::unique(ubait_vec.begin(), ubait_vec.end());
		ubait_vec.resize( tmp_iter - ubait_vec.begin() );
	}
	nbait = ubait_set.size();
	return vector<string>(ubait_set.begin(), ubait_set.end());
}

vector<vector<size_t> > createP2Pmap(const vector<string>& prey_list, const vector<vector<string> >& go_groups) {
	const size_t nprey = prey_list.size();
	map<string, unsigned> prey_index;
	for(unsigned i=0; i<prey_list.size(); ++i){
		prey_index[prey_list[i]] = i;
	}
	Fastmat<unsigned> p2p_matrix(nprey, nprey, 0);
	for(const auto& gene_vector : go_groups)
		for(unsigned i=0; i<gene_vector.size(); ++i)
			for(unsigned j=i+1; j<gene_vector.size(); ++j)
				try {
					unsigned idx1 = prey_index.at(gene_vector[i]);
					unsigned idx2 = prey_index.at(gene_vector[j]);
					++p2p_matrix(idx1, idx2);
					++p2p_matrix(idx2, idx1);
				} catch (std::out_of_range&){}
	vector<vector<size_t> > p2p_mapping(nprey);
	for(size_t r = 0; r < nprey; ++r)
		for(size_t c = 0; c < nprey; ++c)
			if( r != c && p2p_matrix(r, c) >= 1 )
				p2p_mapping[r].push_back(c);
	return p2p_mapping;
}

} // namespace saint_int

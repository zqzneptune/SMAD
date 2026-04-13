/*
 * int_entry.cpp
 * Rcpp entry point for SAINTexpress-int (intensity mode)
 * Replaces original main.cpp + main.hpp file I/O with R DataFrame interface
 */

#include <Rcpp.h>
#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <deque>
#include <iomanip>

// [[Rcpp::depends(BH)]]
#include <boost/array.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/foreach.hpp>

#include "IntGlobals.hpp"
#include "IntStats_fns.hpp"

namespace saint_int {

using namespace std;

// Forward declarations
vector<string> uniqueBait(const deque<BaitClass> &BDATA, size_t &nbait);
vector<vector<size_t> > createP2Pmap(const vector<string>& prey_list, const vector<vector<string> >& go_groups);
Model_data statModel(const Options& opts, const vector<vector<size_t> >& p2p_mapping, const vector<string>& ubait, const Quantmat& test_mat_DATA, const Quantmat& ctrl_mat_DATA, const vector<size_t>& ip_idx_to_bait_no, const size_t nprey, const size_t nbait);

// --- Data parsing from R DataFrames ---

void parseInterFromDF(deque<InterClass>& IDATA, const Rcpp::DataFrame& inter_df, size_t& ninter) {
	Rcpp::StringVector ip_ids = inter_df["ip_id"];
	Rcpp::StringVector bait_ids = inter_df["bait_id"];
	Rcpp::StringVector prey_ids = inter_df["prey_id"];
	Rcpp::NumericVector quants = inter_df["quant"];
	ninter = 0;
	for(int i = 0; i < ip_ids.size(); i++) {
		IDATA.push_back(InterClass());
		InterClass& tmp = IDATA.back();
		tmp.set_ipId(Rcpp::as<string>(ip_ids[i]));
		tmp.set_baitId(Rcpp::as<string>(bait_ids[i]));
		tmp.set_preyId(Rcpp::as<string>(prey_ids[i]));
		// INT mode: log-transform the quantitative value
		tmp.set_quant(log(quants[i]));
		ninter++;
	}
}

set<string> parsePreyFromDF(deque<PreyClass>& PDATA, const Rcpp::DataFrame& prey_df, size_t& nprey) {
	Rcpp::StringVector prey_ids = prey_df["prey_id"];
	bool has_length = prey_df.containsElementNamed("prey_length");
	bool has_gene = prey_df.containsElementNamed("gene_name");
	Rcpp::NumericVector prey_lengths;
	Rcpp::StringVector gene_names;
	if(has_length) prey_lengths = prey_df["prey_length"];
	if(has_gene) gene_names = prey_df["gene_name"];

	nprey = 0;
	set<string> prey_Id_set;
	for(int i = 0; i < prey_ids.size(); i++) {
		string prey_id = Rcpp::as<string>(prey_ids[i]);
		double prey_length = has_length ? prey_lengths[i] : 0.0;
		string gene_name = has_gene ? Rcpp::as<string>(gene_names[i]) : prey_id;

		PreyClass tmp;
		tmp.set_rowId(nprey);
		tmp.set_preyId(prey_id);
		prey_Id_set.insert(prey_id);
		tmp.set_preyLength(prey_length);
		tmp.set_preyGeneId(gene_name);
		PDATA.push_back(tmp);
		nprey++;
	}
	if(nprey != prey_Id_set.size())
		throw runtime_error("duplicate preys in prey data");
	return prey_Id_set;
}

map<string, BaitClass> parseBaitFromDF(deque<BaitClass>& BDATA, const Rcpp::DataFrame& bait_df, size_t& nip, size_t& n_test_ip) {
	Rcpp::StringVector ip_ids = bait_df["ip_id"];
	Rcpp::StringVector bait_ids = bait_df["bait_id"];
	Rcpp::StringVector test_ctrl = bait_df["test_ctrl"];

	nip = 0;
	map<string, BaitClass> bait_Id_map;
	n_test_ip = 0;
	size_t n_ctrl_ip = 0;

	for(int i = 0; i < ip_ids.size(); i++) {
		string node1 = Rcpp::as<string>(ip_ids[i]);
		string node2 = Rcpp::as<string>(bait_ids[i]);
		string node3 = Rcpp::as<string>(test_ctrl[i]);

		BaitClass tmp;
		if(node3 != "T" && node3 != "C")
			throw runtime_error("test_ctrl column must be 'T' or 'C'");
		tmp.set_colId(node3 == "T" ? n_test_ip++ : n_ctrl_ip++);
		tmp.set_ipId(node1);
		tmp.set_baitId(node2);
		tmp.set_isCtrl(node3 != "T");

		BDATA.push_back(tmp);
		bait_Id_map[node1] = tmp;
		nip++;
	}
	return bait_Id_map;
}

void mapRowCol(deque<InterClass>& IDATA, const deque<PreyClass>& PDATA, const deque<BaitClass>& BDATA, const map<string, BaitClass>& bait_Id_map) {
	for(deque<InterClass>::iterator m = IDATA.begin(); m != IDATA.end(); m++) {
		deque<PreyClass>::const_iterator mp = PDATA.begin();
		for(; mp != PDATA.end(); mp++)
			if(mp->get_preyId() == m->get_preyId())
				break;
		if(mp == PDATA.end())
			throw runtime_error("prey " + m->get_preyId() + " not found");
		m->set_rowId(mp->get_rowId());

		deque<BaitClass>::const_iterator mb = BDATA.begin();
		for(; mb != BDATA.end(); mb++)
			if(mb->get_ipId() == m->get_ipId())
				break;
		if(mb == BDATA.end())
			throw runtime_error("bait not found");
		m->set_colId(mb->get_colId());
		m->is_ctrl = bait_Id_map.at(m->get_ipId()).get_isCtrl();
	}
}

int get_nctrl(const deque<BaitClass>& BDATA) {
	int ret = 0;
	for(auto m = BDATA.begin(); m != BDATA.end(); m++)
		if(m->get_isCtrl() == true) ret++;
	return ret;
}

void createMatrixData(Quantmat& test_mat_DATA, Quantmat& ctrl_mat_DATA, const deque<InterClass>& IDATA) {
	for(deque<InterClass>::const_iterator m = IDATA.begin(); m != IDATA.end(); m++) {
		int r = m->get_rowId(), c = m->get_colId();
		Quant_t q = m->get_quant();
		m->is_ctrl ? ctrl_mat_DATA(r, c) = q : test_mat_DATA(r, c) = q;
	}
}

void createList(deque<UIClass>& UIDATA, const deque<InterClass>& IDATA, const deque<BaitClass>& BDATA, const deque<PreyClass>& PDATA, const size_t nprey, const size_t nbait, const vector<string>& ubait, vector<size_t>& ip_idx_to_bait_no, size_t& nuinter) {
	ip_idx_to_bait_no.clear();
	{
		map<string, size_t> ubait_map;
		for(size_t i = 0; i < ubait.size(); ++i)
			ubait_map[ubait[i]] = i;
		BOOST_FOREACH(const auto& bait , BDATA){
			if(!bait.get_isCtrl())
				ip_idx_to_bait_no.push_back(ubait_map.at(bait.get_baitId()));
		}
	}
	UIDATA.clear();
	Fastmat<UIClass*> UImat(nprey, nbait);
	BOOST_FOREACH(const auto& inter , IDATA)
		if(!inter.is_ctrl){
			const size_t j = inter.get_colId();
			const size_t b = ip_idx_to_bait_no[j];
			const size_t i = inter.get_rowId();
			if(!UImat(i,b)){
				UIDATA.push_back({});
				UIClass& tmp = UIDATA.back();
				tmp.set_baitId(inter.get_baitId());
				tmp.set_preyId(inter.get_preyId());
				tmp.set_preyGeneId(PDATA[i].get_preyGeneId());
				tmp.set_rowId(i);
				UImat(i,b) = &tmp;
			}
			UImat(i,b)->add_colId(j);
		}
	nuinter = UIDATA.size();
}

void zero_self_interactions(Quantmat& test_mat_DATA, const deque<UIClass>& UIDATA) {
	BOOST_FOREACH(const auto& UI , UIDATA) {
		if(UI.get_preyGeneId() == UI.get_baitId()) {
			BOOST_FOREACH(const int col , UI.get_colId())
				test_mat_DATA(UI.get_rowId(), col) = saint::nan;
		}
	}
}

// ICM
Model_data icms(Model_data& dp, const bool with_gamma) {
	double oldllik, newllik = dp.llikelihood();
	if(!with_gamma) dp.gamma = 0;else dp.gamma=0;
	for (unsigned iter = 1;iter <= 15; ++iter) {
		dp.print_MRF_parameters();
		oldllik = newllik;
		dp.icm_Z();
		if(with_gamma)
			dp.wrt_MRF();
		else
			dp.wrt_MRF_gamma_0();
		newllik = dp.llikelihood();
		if( (newllik >= oldllik) &&
			(exp(newllik - oldllik) -1 < 1e-3) ) break;
	}
	dp.print_MRF_parameters();
	return dp;
}

boost::array<Fastmat<double>, 3> computation(const bool with_gamma, Model_data& dp, const Options& opts) {
	const auto dp0 = icms(dp, with_gamma);
	const auto all_scores = dp0.calculateScore(opts);
	dp0.print_MRF_parameters();
	return all_scores;
}

} // namespace saint_int


// [[Rcpp::export(.SAINTexpress_int_impl)]]
Rcpp::DataFrame SAINTexpress_int_impl(Rcpp::DataFrame inter_df,
                                       Rcpp::DataFrame prey_df,
                                       Rcpp::DataFrame bait_df,
                                       Rcpp::Nullable<Rcpp::DataFrame> GO_df,
                                       double f, int R, int L) {
	using namespace saint_int;
	using namespace std;

	Options opts;
	opts.f = f;
	opts.R = static_cast<unsigned>(R);
	opts.L = L;

	if(opts.f < 0.0 || 1.0 < opts.f)
		Rcpp::stop("frequency parameter must be between 0 and 1.");
	if(opts.R < 1)
		Rcpp::stop("R must be a positive integer.");
	if(opts.L < 1)
		Rcpp::stop("L must be a positive integer.");

	size_t ninter, nuinter, nbait, nprey, nip;

	// Parse data from R DataFrames
	deque<PreyClass> PDATA;
	parsePreyFromDF(PDATA, prey_df, nprey);

	deque<BaitClass> BDATA;
	size_t n_test_ip;
	const map<string, BaitClass> bait_Id_map = parseBaitFromDF(BDATA, bait_df, nip, n_test_ip);

	deque<InterClass> IDATA;
	parseInterFromDF(IDATA, inter_df, ninter);

	mapRowCol(IDATA, PDATA, BDATA, bait_Id_map);

	// Creating data structure
	const size_t n_ctrl_ip = get_nctrl(BDATA);
	Quantmat test_mat_DATA(nprey, nip - n_ctrl_ip, saint::nan);
	Quantmat ctrl_mat_DATA(nprey, n_ctrl_ip, saint::nan);
	createMatrixData(test_mat_DATA, ctrl_mat_DATA, IDATA);

	// Z-score normalization (INT-specific)
	vector<double> tmp;
	for(auto i=begin(test_mat_DATA.mat); i!=end(test_mat_DATA.mat);++i){
		const Quant_t v = *i;
		if(!saint::isnan(v) && -1e5<v)
			tmp.push_back(v);
	}
	for(auto i=begin(ctrl_mat_DATA.mat); i!=end(ctrl_mat_DATA.mat); ++i){
		const Quant_t v = *i;
		if(!saint::isnan(v) && -1e5<v)
			tmp.push_back(v);
	}
	double tmp_mean = 0;
	for(size_t i = 0; i < tmp.size(); i++) tmp_mean += tmp[i];
	tmp_mean /= tmp.size();
	double tmp_var = 0;
	for(size_t i = 0; i < tmp.size(); i++) tmp_var += (tmp[i] - tmp_mean) * (tmp[i] - tmp_mean);
	tmp_var /= (tmp.size() - 1.0);
	for(size_t i = 0; i < test_mat_DATA.mat.size(); ++i) {
		if(!saint::isnan(test_mat_DATA.mat[i])) {
			test_mat_DATA.mat[i] = (static_cast<double>(test_mat_DATA.mat[i]) - tmp_mean) / sqrt(tmp_var);
		}
	}
	for(size_t i = 0; i < ctrl_mat_DATA.mat.size(); ++i) {
		if(!saint::isnan(ctrl_mat_DATA.mat[i])) {
			ctrl_mat_DATA.mat[i] = (static_cast<double>(ctrl_mat_DATA.mat[i]) - tmp_mean) / sqrt(tmp_var);
		}
	}

	// Unique interactions
	const vector<string> ubait = uniqueBait(BDATA, nbait);
	vector<size_t> ip_idx_to_bait_no;
	deque<UIClass> UIDATA;
	createList(UIDATA, IDATA, BDATA, PDATA, nprey, nbait, ubait, ip_idx_to_bait_no, nuinter);
	zero_self_interactions(test_mat_DATA, UIDATA);

	// Mapping proteins
	vector<string> prey_list;
	for(size_t i = 0; i<nprey; ++i)
		prey_list.push_back(PDATA.at(i).get_preyId());

	// Parse GO data if provided
	vector<vector<size_t> > p2p_mapping;
	bool has_GO = GO_df.isNotNull();
	if(has_GO) {
		Rcpp::DataFrame go = Rcpp::as<Rcpp::DataFrame>(GO_df);
		Rcpp::StringVector go_genes_col = go["genes"];
		vector<vector<string> > go_groups;
		for(int i = 0; i < go_genes_col.size(); i++) {
			string genes_str = Rcpp::as<string>(go_genes_col[i]);
			vector<string> gene_vector;
			boost::algorithm::split(gene_vector, genes_str, boost::is_any_of(" "));
			go_groups.push_back(gene_vector);
		}
		p2p_mapping = createP2Pmap(prey_list, go_groups);
	} else {
		p2p_mapping = vector<vector<size_t> >(prey_list.size(), vector<size_t>(0));
	}

	// Statistical analysis
	Model_data dp = statModel(opts, p2p_mapping, ubait, test_mat_DATA, ctrl_mat_DATA, ip_idx_to_bait_no, nprey, nbait);

	Model_data dp_gamma0 = dp;
	const boost::array<Fastmat<double>, 3> scores = computation(false, dp_gamma0, opts);
	Model_data dp_MRF = dp;
	const auto ctor = Fastmat<double>(nprey, nbait);
	boost::array<Fastmat<double>, 3> topo_scores{{ctor, ctor, ctor}};
	if(has_GO)
		topo_scores = computation(true, dp_MRF, opts);
	else
		topo_scores = scores;

	// Build output DataFrame
	Rcpp::StringVector out_Bait, out_Prey, out_PreyGene, out_Intensity, out_ctrlIntensity, out_boosted_by;
	Rcpp::NumericVector out_IntensitySum, out_AvgIntensity, out_NumReplicates;
	Rcpp::NumericVector out_AvgP, out_MaxP, out_TopoAvgP, out_TopoMaxP;
	Rcpp::NumericVector out_SaintScore, out_OddsScore, out_FoldChange, out_BFDR;

	BOOST_FOREACH(const auto& ui, UIDATA) {
		size_t row = ui.get_rowId();
		size_t baitcol = ip_idx_to_bait_no[ui.get_colId().front()];

		out_Bait.push_back(ui.get_baitId());
		out_Prey.push_back(ui.get_preyId());
		out_PreyGene.push_back(ui.get_preyGeneId());

		// Original intensities (exp transform back from log)
		valarray<Quant_t> test_row = dp_MRF.test_mat(row, baitcol);
		valarray<Quant_t> quant_orig(test_row.size());
		for(size_t k=0; k<test_row.size(); ++k) quant_orig[k] = exp(test_row[k]);
		valarray<Quant_t> quant = quant_orig;
		for(size_t i=0; i<quant.size(); ++i)
			if(saint::isnan(quant[i]))
				quant[i] = exp(dp_MRF.t);
		double AvgQuant = quant.sum()/quant.size();

		// Build pipe-separated intensity string (with NaN as ".")
		string int_str;
		{
			std::ostringstream oss;
			oss << std::setprecision(3) << std::fixed;
			if(saint::isnan(quant_orig[0])) oss << ".";
			else oss << static_cast<double>(quant_orig[0]);
			for(size_t col = 1; col < quant_orig.size(); ++col) {
				oss << "|";
				if(saint::isnan(quant_orig[col])) oss << ".";
				else oss << static_cast<double>(quant_orig[col]);
			}
			int_str = oss.str();
		}
		out_Intensity.push_back(int_str);

		out_IntensitySum.push_back(quant.sum());
		out_AvgIntensity.push_back(AvgQuant);
		out_NumReplicates.push_back(quant.size());

		// Control intensities
		valarray<Quant_t> ctrl_row = ctrl_mat_DATA(row, __);
		valarray<Quant_t> ctrlCounts(ctrl_row.size());
		for(size_t k=0; k<ctrl_row.size(); ++k) ctrlCounts[k] = exp(ctrl_row[k]);
		string ctrl_str;
		{
			std::ostringstream oss;
			oss << std::setprecision(3) << std::fixed;
			if(saint::isnan(ctrlCounts[0])) oss << ".";
			else oss << static_cast<double>(ctrlCounts[0]);
			for(size_t j = 1; j < ctrlCounts.size(); ++j) {
				oss << "|";
				if(saint::isnan(ctrlCounts[j])) oss << ".";
				else oss << static_cast<double>(ctrlCounts[j]);
			}
			ctrl_str = oss.str();
		}
		out_ctrlIntensity.push_back(ctrl_str);

		double avg_p = scores[0](row, baitcol);
		unsigned denom = 0;
		double numer = 0;
		for(unsigned i = 0; i < scores[0].size(); ++i) {
			bool tmp_b = scores[0][i] > avg_p;
			numer += scores[0][i] * tmp_b;
			denom += tmp_b;
		}
		double BFDR = (denom == 0 ? 0 : 1 - numer / denom);
		double max_p = scores[1](row, baitcol);
		double topo_avg_p = topo_scores[0](row, baitcol);
		double topo_max_p = topo_scores[1](row, baitcol);

		out_AvgP.push_back(avg_p);
		out_MaxP.push_back(max_p);
		out_TopoAvgP.push_back(topo_avg_p);
		out_TopoMaxP.push_back(topo_max_p);
		out_SaintScore.push_back(max(topo_avg_p, avg_p));
		out_OddsScore.push_back(topo_scores[2](row, baitcol));
		out_FoldChange.push_back(AvgQuant / std::exp(dp_MRF.eta[row]));
		out_BFDR.push_back(BFDR);

		// Boosted by
		string boosted = "";
		BOOST_FOREACH(const unsigned l, p2p_mapping[row]) {
			if(dp_MRF.Z(l, baitcol).mean() > 0)
				boosted += prey_list[l] + "|";
		}
		out_boosted_by.push_back(boosted);
	}

	return Rcpp::DataFrame::create(
		Rcpp::Named("Bait") = out_Bait,
		Rcpp::Named("Prey") = out_Prey,
		Rcpp::Named("PreyGene") = out_PreyGene,
		Rcpp::Named("Intensity") = out_Intensity,
		Rcpp::Named("IntensitySum") = out_IntensitySum,
		Rcpp::Named("AvgIntensity") = out_AvgIntensity,
		Rcpp::Named("NumReplicates") = out_NumReplicates,
		Rcpp::Named("ctrlIntensity") = out_ctrlIntensity,
		Rcpp::Named("AvgP") = out_AvgP,
		Rcpp::Named("MaxP") = out_MaxP,
		Rcpp::Named("TopoAvgP") = out_TopoAvgP,
		Rcpp::Named("TopoMaxP") = out_TopoMaxP,
		Rcpp::Named("SaintScore") = out_SaintScore,
		Rcpp::Named("OddsScore") = out_OddsScore,
		Rcpp::Named("FoldChange") = out_FoldChange,
		Rcpp::Named("BFDR") = out_BFDR,
		Rcpp::Named("boosted_by") = out_boosted_by,
		Rcpp::Named("stringsAsFactors") = false
	);
}

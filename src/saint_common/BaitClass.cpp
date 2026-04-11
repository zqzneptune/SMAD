#include <Rcpp.h>
/*
 * bait_class.cpp
 *
 *  Created on: Jun 18, 2012
 *      Author: hwchoi
 */

#include <iostream>
#include <string>
#include <fstream>
#include <vector>

#include "BaitClass.hpp"

using namespace std;


void BaitClass::print() const {

	Rcpp::Rcout << colId << "\t";
	Rcpp::Rcout << ipId << "\t";
	Rcpp::Rcout << baitId << "\t";
	Rcpp::Rcout << isCtrl << endl;

}

int BaitClass::get_colId() const {
	return colId;
}

string BaitClass::get_ipId() const {
	return ipId;
}

const string & BaitClass::get_baitId() const {
	return baitId;
}

bool BaitClass::get_isCtrl() const {
	return isCtrl;
}

void BaitClass::set_colId( int c ) {
	colId = c;
}

void BaitClass::set_ipId( string str ) {
	ipId = str;
}

void BaitClass::set_baitId( string str ) {
	baitId = str;
}

void BaitClass::set_isCtrl( bool ic ) {
	isCtrl = ic;
}

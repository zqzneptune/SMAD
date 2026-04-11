#include <Rcpp.h>
/*
 * prey_class.cpp
 *
 *  Created on: Jun 18, 2012
 *      Author: hwchoi
 */

#include <iostream>
#include <string>
#include <fstream>
#include <vector>

#include "PreyClass.hpp"

using namespace std;

PreyClass::PreyClass() {

}


void PreyClass::print() const {

	Rcpp::Rcout << rowId << "\t";
	Rcpp::Rcout << preyId << "\t";
	Rcpp::Rcout << (int) preyLength << "\t";
	Rcpp::Rcout << preyGeneId << endl;

}


int PreyClass::get_rowId() const {
	return rowId;
}

string PreyClass::get_preyId() const {
	return preyId;
}

double PreyClass::get_preyLength() const {
	return preyLength;
}

string PreyClass::get_preyGeneId() const {
	return preyGeneId;
}

void PreyClass::set_rowId(int r) {
	rowId = r;
}

void PreyClass::set_preyId(string str) {
	preyId = str;
}

void PreyClass::set_preyLength(double x) {
	preyLength = x;
}

void PreyClass::set_preyGeneId(string str) {
	preyGeneId = str;
}

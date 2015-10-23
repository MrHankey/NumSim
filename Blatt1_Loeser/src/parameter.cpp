#include "parameter.hpp"

#include <iostream>
#include <sstream>
#include <fstream>
using namespace std;

Parameter::Parameter() {
	_re = 1000.0;
	_omega = 1.7;
	_alpha = 0.9;
	_dt = 0.2;
	_tend = 16.4;
	_eps = 0.001;
	_tau = 0.5;
	_itermax = 100;
}

void Parameter::Load(const char* file) {

	ifstream fileStream(file);

	if ( fileStream.fail() )
	{
		cout << "Couldn't load " << file << ". Using default values." << endl;
	}
	else
	{
		fileStream >> _re;
		fileStream >> _omega;
		fileStream >> _alpha;
		fileStream >> _dt;
		fileStream >> _tend;
		fileStream >> _eps;
		fileStream >> _tau;
		fileStream >> _itermax;
	}

	cout << "Loaded parameters from " << file << ":" << endl;
	cout << _re << endl;
	cout << _omega << endl;
	cout << _alpha << endl;
	cout << _dt << endl;
	cout << _tend << endl;
	cout << _eps << endl;
	cout << _tau << endl;
	cout << _itermax << endl;

}

const real_t& Parameter::Re() const {
	return _re;
}

const real_t& Parameter::Omega() const {
	return _omega;
}

const real_t& Parameter::Alpha() const {
	return _alpha;
}

const real_t& Parameter::Dt() const {
	return _dt;
}

const real_t& Parameter::Tend() const {
	return _tend;
}

const index_t& Parameter::IterMax() const {
	return _itermax;
}

const real_t& Parameter::Eps() const {
	return _eps;
}

const real_t& Parameter::Tau() const {
	return _tau;
}

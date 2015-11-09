/*
 * Copyright (C) 2015  Raphael Leiteriz, Sebastian Reuschen, Hamzeh Kraus
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "parameter.hpp"

#include <iostream>
#include <sstream>
#include <fstream>
using namespace std;

/// Constructs a new Parameter set with default values
Parameter::Parameter() {
	_re      = 1000.0;
	_omega   = 1.7;
	_alpha   = 0.9;
	_dt      = 0.2;
	_tend    = 16.4;
	_eps     = 0.001;
	_tau     = 0.5;
	_itermax = 100;

	cout << "Loaded default parameters"<< endl;
	PrintVariables();
}

/// Loads the parameter values from a file
//  @param file  load data from default file
void Parameter::Load(const char* file) {
	// Initialize
	ifstream fileStream(file);

	// Catch lodaing error
	if ( fileStream.fail() ) {
		cout << "Couldn't load parameters from " << file << endl;
	} else {
		fileStream >> _re;
		fileStream >> _omega;
		fileStream >> _alpha;
		fileStream >> _dt;
		fileStream >> _tend;
		fileStream >> _eps;
		fileStream >> _tau;
		fileStream >> _itermax;
		cout << "Loaded parameters from " << file << ":" << endl;
		PrintVariables();
	}
}

///Prints Parameters
void Parameter::PrintVariables(){
	cout << "Re: "      << _re      << endl;
	cout << "omega: "   << _omega   << endl;
	cout << "alpha: "   << _alpha   << endl;
	cout << "dt: "      << _dt      << endl;
	cout << "Tend: "    << _tend    << endl;
	cout << "eps: "     << _eps     << endl;
	cout << "tau: "     << _tau     << endl;
	cout << "itermax: " << _itermax << endl << endl;
}


/* Getter functions for all parameters */
const real_t&  Parameter::Re()      const {return _re;}
const real_t&  Parameter::Omega()   const {return _omega;}
const real_t&  Parameter::Alpha()   const {return _alpha;}
const real_t&  Parameter::Dt()      const {return _dt;}
const real_t&  Parameter::Tend()    const {return _tend;}
const index_t& Parameter::IterMax() const {return _itermax;}
const real_t&  Parameter::Eps()     const {return _eps;}
const real_t&  Parameter::Tau()     const {return _tau;}

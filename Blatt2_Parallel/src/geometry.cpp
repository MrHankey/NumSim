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

#include "geometry.hpp"
#include "grid.hpp"
#include "iterator.hpp"
#include "communicator.hpp"

#include <fstream>
#include <iostream>
#include <stdexcept>
#include <sstream>

using namespace std;

/// Constructs a default geometry:
/// driven cavity with 128 x 128 grid, no-slip boundary conditions
Geometry::Geometry() {
	// Set starting velocities
	_velocity[0] = 1;
	_velocity[1] = 0;
	_pressure    = 0;

	_comm = nullptr;

	// Number of cells in one line
	_bsize[0] = 128;
	_bsize[1] = 128;

	// Length of driven cavity
	_blength[0] = 1;
	_blength[1] = 1;

	// Print vars
	cout << "Loaded default geometry definition." << endl;
}

Geometry::Geometry(Communicator* comm) : Geometry() {
	_comm = comm;
	Initialize();
}

/// Loads a geometry from a file
//  @param file  txt file to load data from
void Geometry::Load(const char* file) {
	// Initialize
	ifstream fileStream(file);

	// Cath loading error
	if ( fileStream.fail() ) {
		// Failed to load
		cout << "Couldn't load geometry from " << file << endl;
		return;
	} else {
		//Load values
		fileStream >> _velocity[0];
		fileStream >> _velocity[1];
		fileStream >> _pressure;
		fileStream >> _bsize[0];
		fileStream >> _bsize[1];
		fileStream >> _blength[0];
		fileStream >> _blength[1];

		// Success
		cout << "Loaded geometry definitions from " << file << "." << endl;
	}
	Initialize();
}

void Geometry::Initialize()
{
	// Cell Length
	_h[0] = _blength[0]/(_bsize[0]);
	_h[1] = _blength[1]/(_bsize[1]);

	_size[0] = _bsize[0]/_comm->ThreadDim()[0] + 2;
	_size[1] = _bsize[1]/_comm->ThreadDim()[1] + 2;

	if ( _bsize[0] % _comm->ThreadDim()[0] != 0 || _bsize[1] % _comm->ThreadDim()[1] != 0 )
	{
		throw std::runtime_error("Number of cells not devisible by processor distribution");
	}

	_length[0] = _blength[0]/_comm->ThreadDim()[0];
	_length[1] = _blength[1]/_comm->ThreadDim()[1];

	_bsize[0] += 2;
	_bsize[1] += 2;

	if ( (_bsize[0] % 2) == 0 && (_bsize[1] % 2) == 0)
	{
		_comm->SetEvenOdd(true);
	}
	else if (_bsize[0] % 2 != 0 && _bsize[1] % 2 != 0)
	{
		_comm->SetEvenOdd((_comm->ThreadIdx()[0] + _comm->ThreadIdx()[1]) % 2 == 0 );
	}
	else if (_bsize[0] % 2 != 0 && _bsize[1] % 2 == 0)
	{
		_comm->SetEvenOdd(_comm->ThreadIdx()[0] % 2 == 0);
	}
	else if (_bsize[0] % 2 == 0 && _bsize[1] % 2 != 0)
	{
		_comm->SetEvenOdd(_comm->ThreadIdx()[1] % 2 == 0);
	}

	/*multi_index_t procPos = _comm->ThreadIdx();
	multi_index_t startCell = {procPos[0]*(_size[0]) , procPos[1]*(_size[1])};
	bool evenodd = (startCell[0] + _bsize[0]*startCell[1]) != 0;
	_comm->SetEvenOdd(evenodd);*/

	//read CSV file
	_b  = new Grid(this);
	this->readCsvGrid("yolo.swag");

	printf(" local_siz: %i \n", _size[0]);
}

/// Prints Parameters
void Geometry::PrintVariables(){
	cout << "vel_x: "    << _velocity[0]  << endl;
	cout << "vel_y: "    << _velocity[1]  << endl;
	cout << "p: "        << _pressure     << endl;
	cout << "size_x: "   << _bsize[0]      << endl;
	cout << "size_y: "   << _bsize[1]      << endl;
	cout << "length_x: " << _blength[0]    << endl;
	cout << "length_y "  << _blength[1]    << endl;
	cout << "h_x: "      << _h[0]         << endl;
	cout << "h_y: "      << _h[1] << endl << endl;

}

const multi_index_t& Geometry::TotalSize() const {
	return _bsize;
}

const multi_real_t& Geometry::TotalLength() const {
	return _blength;
}

/* Getter functions */
/// Returns the number of cells in each dimension
const multi_index_t& Geometry::Size()   const {
	return _size;
}
/// Returns the length of the domain
const multi_real_t&  Geometry::Length() const
{
	return _length;
}
/// Returns the meshwidth
const multi_real_t&  Geometry::Mesh()   const {return _h;}


/// Read Csv and set in boundary grid
// @param fileName  name of file to read
void Geometry::readCsvGrid(string fileName) const {
	ifstream file(fileName);
	string value;
	Iterator it = Iterator(this);
	string line;
	while(getline(file, line)) {

		while ( line.find ("\r\n") != string::npos )
		{
			line.erase ( line.find ("\r\n"), 2 );
		}

		while ( value.find ("\n") != string::npos )
		{
			value.erase ( value.find ("\n"), 2 );
		}

		std::stringstream lineStream(line);

		while ( getline(lineStream, value, ',') )
		{
			while ( line.find ("\r\n") != string::npos )
			{
				line.erase ( line.find ("\r\n"), 2 );
			}

			while ( value.find ("\n") != string::npos )
			{
				value.erase ( value.find ("\n"), 2 );
			}


			_b->Cell(it) = atoi(value.c_str());
			cout << _b->Cell(it);
			//cout << atoi(value.c_str());
			it.Next();
		}
	}

	cout<<endl;
	Grid *_bTmp;
	_bTmp = new Grid(this);
	Iterator ita = Iterator(this);
	while(ita.Valid()) {
		_bTmp->Cell(ita) = _b->Cell(ita);
		ita.Next();
	}

	Iterator its = Iterator(this);
	while(its.Valid()) {
		index_t numCol  = its.Value()%this->_bsize[0];
		index_t numRow  = (its.Value())/this->_bsize[0];
		index_t flipRow = (index_t)abs((int)numRow-(int)this->_bsize[1])-1;
		index_t newCell = flipRow*this->_bsize[0]+numCol;
		Iterator newPos = Iterator(this,newCell);

		//cout << newPos.Value() <<", " << newCell <<"; ";

		cout << _bTmp->Cell(its);

		_b->Cell(newPos)=_bTmp->Cell(its);
		//cout << _b->Cell(newPos);
		its.Next();
	}
}


/// Updates the boundary velocity field u
//  @param u  get grid of u
void Geometry::Update_U(Grid* u) const {
	// Initialize
	real_t velocity = _velocity[0];
	BoundaryIterator it = BoundaryIterator(this);

	// Bottom boundary update
	if(_comm->isBottom()){
		it.SetBoundary(it.boundaryBottom);
		while ( it.Valid() ) {
			u->Cell(it) = -1*u->Cell(it.Top());
			it.Next();
		}
	}

	// Right boundary update
	if(_comm->isRight()){
		it.SetBoundary(it.boundaryRight);
		while ( it.Valid() ) {
			u->Cell(it) = 0;
			u->Cell(it.Left()) = 0; // hier evtl deleten
			it.Next();
		}
	}
	// Top boundary update
	if(_comm->isTop()){
		it.SetBoundary(it.boundaryTop);
		while ( it.Valid() ) {
			u->Cell(it) = 2*velocity - u->Cell(it.Down()) ;
			it.Next();
		}
	}
	// Left boundary update
	if(_comm->isLeft()){
		it.SetBoundary(it.boundaryLeft);
		while ( it.Valid() ) {
			u->Cell(it) = 0;
			it.Next();
		}
	}
}

/// Updates the boundary velocity field v
//  @param v  get grid of v
void Geometry::Update_V(Grid* v) const {
	// Initialize
	BoundaryIterator it = BoundaryIterator(this);

	// Bottom boundary update
	if(_comm->isBottom()){
		it.SetBoundary(it.boundaryBottom);
		while ( it.Valid() ) {
			v->Cell(it) = 0;
			it.Next();
		}
	}
	// Right boundary update
	if(_comm->isRight()){
		it.SetBoundary(it.boundaryRight);
		while ( it.Valid() ) {
			v->Cell(it) = -1*v->Cell(it.Left());
			it.Next();
		}
	}
	// Top boundary update
	if(_comm->isTop()){
		it.SetBoundary(it.boundaryTop);
		while ( it.Valid() ) {
			v->Cell(it.Down()) = 0; // hier evtl deleten
			v->Cell(it) = 0;
			it.Next();
		}
	}
	// Left boundary update
	if(_comm->isLeft()){
		it.SetBoundary(it.boundaryLeft);
		while ( it.Valid() ) {
			v->Cell(it) = -1*v->Cell(it.Right());
			it.Next();
		}
	}
}

/// Updates the boundary velocity field p
//  @param p  get grid of p
void Geometry::Update_P(Grid* p) const {
	// Initialize
	BoundaryIterator it = BoundaryIterator(this);
	it.SetBoundary(it.boundaryBottom);

	// Bottom boundary update
	if(_comm->isBottom()){
		while ( it.Valid() ) {
			p->Cell(it) = p->Cell(it.Top());
			it.Next();
		}
	}
	// Right boundary update
	if(_comm->isRight()){
		it.SetBoundary(it.boundaryRight);
		while ( it.Valid() ) {
			p->Cell(it) = p->Cell(it.Left());
			it.Next();
		}
	}
	// Top boundary update
	if(_comm->isTop()){
		it.SetBoundary(it.boundaryTop);
		while ( it.Valid() ) {
			p->Cell(it) = p->Cell(it.Down());
			it.Next();
		}
	}
	// Left boundary update
	if(_comm->isLeft()){
		it.SetBoundary(it.boundaryLeft);
		while ( it.Valid() ) {
			p->Cell(it) = p->Cell(it.Right());
			it.Next();
		}
	}
}

void Geometry::Update_All(Grid* p,Grid* u,Grid* v,real_t pL,real_t pR) const {
	Iterator it = Iterator(this);
	it.First();
	while (it.Valid()){
		//cout<<it.Value()<<"b:"<< _b->Cell(it)<<endl;
		if(_b->Cell(it)==0){
			//Mach nichts da keine Randpunkt
		}
		else if(_b->Cell(it)==1){//no Slip Bedingung
			if(_b->Cell(it.Right())!=0 && _b->Cell(it.Top())!=0 && _b->Cell(it.Left())!=0 && _b->Cell(it.Down())!=0){
				//nichts tun da alle Ränder auch Rand sind.
			}
			else if(_b->Cell(it.Right())==0 && _b->Cell(it.Top())!=0 && _b->Cell(it.Left())!=0 && _b->Cell(it.Down())!=0){
				//Rechts ist Wasser
				p->Cell(it)=p->Cell(it.Right());
				u->Cell(it)=0;
				v->Cell(it)=-1*(v->Cell(it.Right()));
			}
			else if(_b->Cell(it.Right())!=0 && _b->Cell(it.Top())==0 && _b->Cell(it.Left())!=0 && _b->Cell(it.Down())!=0){
				//Oben ist Wasser
				p->Cell(it)=p->Cell(it.Top());
				u->Cell(it)=-1*(u->Cell(it.Top()));
				v->Cell(it)=0;
			}
			else if(_b->Cell(it.Right())!=0 && _b->Cell(it.Top())!=0 && _b->Cell(it.Left())==0 && _b->Cell(it.Down())!=0){
				//Links ist Wasser
				p->Cell(it)=p->Cell(it.Left());
				u->Cell(it)=0;
				u->Cell(it.Left())=0;
				v->Cell(it)=-1*(v->Cell(it.Left()));
			}
			else if(_b->Cell(it.Right())!=0 && _b->Cell(it.Top())!=0 && _b->Cell(it.Left())!=0 && _b->Cell(it.Down())==0){
				//Unten ist Wasser
				p->Cell(it)=p->Cell(it.Down());
				u->Cell(it)=-1*(u->Cell(it.Down()));
				v->Cell(it)=0;
				v->Cell(it.Down())=0;
			}
			else if(_b->Cell(it.Right())==0 && _b->Cell(it.Top())==0 && _b->Cell(it.Left())!=0 && _b->Cell(it.Down())!=0){
				//Rechts und Oben ist Wasser
				p->Cell(it)=( p->Cell(it.Right())+p->Cell(it.Top()) )/2; // Mittelwert des Druckes
				u->Cell(it)=0;
				v->Cell(it)=0;
			}
			else if(_b->Cell(it.Right())!=0 && _b->Cell(it.Top())==0 && _b->Cell(it.Left())==0 && _b->Cell(it.Down())!=0){
				//Oben und Links ist Wasser
				p->Cell(it)=( p->Cell(it.Left())+p->Cell(it.Top()) )/2; // Mittelwert des Druckes
				u->Cell(it)=0;
				u->Cell(it.Left());
				v->Cell(it)=0;
			}
			else if(_b->Cell(it.Right())!=0 && _b->Cell(it.Top())!=0 && _b->Cell(it.Left())==0 && _b->Cell(it.Down())==0){
				//Links und Unten ist Wasser
				p->Cell(it)=( p->Cell(it.Left())+p->Cell(it.Down()) )/2; // Mittelwert des Druckes
				u->Cell(it)=0;
				u->Cell(it.Left());
				v->Cell(it)=0;
				v->Cell(it.Down())=0;

			}
			else if(_b->Cell(it.Right())==0 && _b->Cell(it.Top())!=0 && _b->Cell(it.Left())!=0 && _b->Cell(it.Down())==0){
				//Unten und Rechts ist Wasser
				p->Cell(it)=( p->Cell(it.Down())+p->Cell(it.Right()) )/2; // Mittelwert des Druckes
				u->Cell(it)=0;
				v->Cell(it)=0;
				v->Cell(it.Down())=0;
			}
			else{
				cout<<"Error in CSV File"<< it.Value()<<"Bedingung:"<< _b->Cell(it)<<endl;
			}

		}
		else if(_b->Cell(it)==2){//Slip Bedingung
			cout<<"2 not implemented!"<<endl;
			if(_b->Cell(it.Right())!=0 && _b->Cell(it.Top())!=0 && _b->Cell(it.Left())!=0 && _b->Cell(it.Down())!=0){
				//nichts tun da alle Ränder auch Rand sind.
			}
			else if(_b->Cell(it.Right())==0 && _b->Cell(it.Top())!=0 && _b->Cell(it.Left())!=0 && _b->Cell(it.Down())!=0){
				//Rechts ist Wasser
			}
			else if(_b->Cell(it.Right())!=0 && _b->Cell(it.Top())==0 && _b->Cell(it.Left())!=0 && _b->Cell(it.Down())!=0){
				//Oben ist Wasser
			}
			else if(_b->Cell(it.Right())!=0 && _b->Cell(it.Top())!=0 && _b->Cell(it.Left())==0 && _b->Cell(it.Down())!=0){
				//Links ist Wasser
			}
			else if(_b->Cell(it.Right())!=0 && _b->Cell(it.Top())!=0 && _b->Cell(it.Left())!=0 && _b->Cell(it.Down())==0){
				//Unten ist Wasser
			}
			else if(_b->Cell(it.Right())==0 && _b->Cell(it.Top())==0 && _b->Cell(it.Left())!=0 && _b->Cell(it.Down())!=0){
				//Rechts und Oben ist Wasser
			}
			else if(_b->Cell(it.Right())!=0 && _b->Cell(it.Top())==0 && _b->Cell(it.Left())==0 && _b->Cell(it.Down())!=0){
				//Oben und Links ist Wasser
			}
			else if(_b->Cell(it.Right())!=0 && _b->Cell(it.Top())!=0 && _b->Cell(it.Left())==0 && _b->Cell(it.Down())==0){
				//Links und Unten ist Wasser
			}
			else if(_b->Cell(it.Right())==0 && _b->Cell(it.Top())!=0 && _b->Cell(it.Left())!=0 && _b->Cell(it.Down())==0){
				//Unten und Rechts ist Wasser
			}
			else{
				cout<<"Error in CSV File"<< it.Value()<<"Bedingung:"<< _b->Cell(it)<<endl;
			}


		}
		else if(_b->Cell(it)==3){//linker Rand OUTFLOW Bedingung
			if(_b->Cell(it.Right())!=0 && _b->Cell(it.Top())!=0 && _b->Cell(it.Left())!=0 && _b->Cell(it.Down())!=0){
				//nichts tun da alle Ränder auch Rand sind.
			}
			else if(_b->Cell(it.Right())==0 && _b->Cell(it.Top())!=0 && _b->Cell(it.Left())!=0 && _b->Cell(it.Down())!=0){
				//Rechts ist Wasser
				p->Cell(it)= pL;
				u->Cell(it)=u->Cell(it.Right());
				v->Cell(it)=v->Cell(it.Right());
			}
			else if(_b->Cell(it.Right())!=0 && _b->Cell(it.Top())==0 && _b->Cell(it.Left())!=0 && _b->Cell(it.Down())!=0){
				//Oben ist Wasser
			}
			else if(_b->Cell(it.Right())!=0 && _b->Cell(it.Top())!=0 && _b->Cell(it.Left())==0 && _b->Cell(it.Down())!=0){
				//Links ist Wasser
			}
			else if(_b->Cell(it.Right())!=0 && _b->Cell(it.Top())!=0 && _b->Cell(it.Left())!=0 && _b->Cell(it.Down())==0){
				//Unten ist Wasser
			}
			else if(_b->Cell(it.Right())==0 && _b->Cell(it.Top())==0 && _b->Cell(it.Left())!=0 && _b->Cell(it.Down())!=0){
				//Rechts und Oben ist Wasser
			}
			else if(_b->Cell(it.Right())!=0 && _b->Cell(it.Top())==0 && _b->Cell(it.Left())==0 && _b->Cell(it.Down())!=0){
				//Oben und Links ist Wasser
			}
			else if(_b->Cell(it.Right())!=0 && _b->Cell(it.Top())!=0 && _b->Cell(it.Left())==0 && _b->Cell(it.Down())==0){
				//Links und Unten ist Wasser
			}
			else if(_b->Cell(it.Right())==0 && _b->Cell(it.Top())!=0 && _b->Cell(it.Left())!=0 && _b->Cell(it.Down())==0){
				//Unten und Rechts ist Wasser
			}
			else{
				cout<<"Error in CSV File"<< it.Value()<<"Bedingung:"<< _b->Cell(it)<<endl;
			}


		}
		else if(_b->Cell(it)==4){//rechter Rand OUTFLOW Bedingung
			if(_b->Cell(it.Right())!=0 && _b->Cell(it.Top())!=0 && _b->Cell(it.Left())!=0 && _b->Cell(it.Down())!=0){
				//nichts tun da alle Ränder auch Rand sind.
			}
			else if(_b->Cell(it.Right())==0 && _b->Cell(it.Top())!=0 && _b->Cell(it.Left())!=0 && _b->Cell(it.Down())!=0){
				//Rechts ist Wasser
			}
			else if(_b->Cell(it.Right())!=0 && _b->Cell(it.Top())==0 && _b->Cell(it.Left())!=0 && _b->Cell(it.Down())!=0){
				//Oben ist Wasser
			}
			else if(_b->Cell(it.Right())!=0 && _b->Cell(it.Top())!=0 && _b->Cell(it.Left())==0 && _b->Cell(it.Down())!=0){
				//Links ist Wasser
				p->Cell(it)= pR;
				u->Cell(it)=u->Cell(it.Left());
				v->Cell(it)=v->Cell(it.Left());
			}
			else if(_b->Cell(it.Right())!=0 && _b->Cell(it.Top())!=0 && _b->Cell(it.Left())!=0 && _b->Cell(it.Down())==0){
				//Unten ist Wasser
			}
			else if(_b->Cell(it.Right())==0 && _b->Cell(it.Top())==0 && _b->Cell(it.Left())!=0 && _b->Cell(it.Down())!=0){
				//Rechts und Oben ist Wasser
			}
			else if(_b->Cell(it.Right())!=0 && _b->Cell(it.Top())==0 && _b->Cell(it.Left())==0 && _b->Cell(it.Down())!=0){
				//Oben und Links ist Wasser
			}
			else if(_b->Cell(it.Right())!=0 && _b->Cell(it.Top())!=0 && _b->Cell(it.Left())==0 && _b->Cell(it.Down())==0){
				//Links und Unten ist Wasser
			}
			else if(_b->Cell(it.Right())==0 && _b->Cell(it.Top())!=0 && _b->Cell(it.Left())!=0 && _b->Cell(it.Down())==0){
				//Unten und Rechts ist Wasser
			}
			else{
				cout<<"Error in CSV File"<< it.Value()<<"Bedingung:"<< _b->Cell(it)<<endl;
			}


		}
		else if(_b->Cell(it)==5){//oben INFLOW Bedingung (driven Cavity)
			if(_b->Cell(it.Right())!=0 && _b->Cell(it.Top())!=0 && _b->Cell(it.Left())!=0 && _b->Cell(it.Down())!=0){
				//nichts tun da alle Ränder auch Rand sind.
			}
			else if(_b->Cell(it.Right())==0 && _b->Cell(it.Top())!=0 && _b->Cell(it.Left())!=0 && _b->Cell(it.Down())!=0){
				//Rechts ist Wasser
			}
			else if(_b->Cell(it.Right())!=0 && _b->Cell(it.Top())==0 && _b->Cell(it.Left())!=0 && _b->Cell(it.Down())!=0){
				//Oben ist Wasser
			}
			else if(_b->Cell(it.Right())!=0 && _b->Cell(it.Top())!=0 && _b->Cell(it.Left())==0 && _b->Cell(it.Down())!=0){
				//Links ist Wasser
			}
			else if(_b->Cell(it.Right())!=0 && _b->Cell(it.Top())!=0 && _b->Cell(it.Left())!=0 && _b->Cell(it.Down())==0){
				//Unten ist Wasser
				p->Cell(it)=p->Cell(it.Down());
				u->Cell(it)=2-(u->Cell(it.Down()));
				v->Cell(it)=0;
				v->Cell(it.Down())=0;

			}
			else if(_b->Cell(it.Right())==0 && _b->Cell(it.Top())==0 && _b->Cell(it.Left())!=0 && _b->Cell(it.Down())!=0){
				//Rechts und Oben ist Wasser
			}
			else if(_b->Cell(it.Right())!=0 && _b->Cell(it.Top())==0 && _b->Cell(it.Left())==0 && _b->Cell(it.Down())!=0){
				//Oben und Links ist Wasser
			}
			else if(_b->Cell(it.Right())!=0 && _b->Cell(it.Top())!=0 && _b->Cell(it.Left())==0 && _b->Cell(it.Down())==0){
				//Links und Unten ist Wasser
			}
			else if(_b->Cell(it.Right())==0 && _b->Cell(it.Top())!=0 && _b->Cell(it.Left())!=0 && _b->Cell(it.Down())==0){
				//Unten und Rechts ist Wasser
			}
			else{
				cout<<"Error in CSV File"<< it.Value()<<"Bedingung:"<< _b->Cell(it)<<endl;
			}


		}
		it.Next();
	}
}

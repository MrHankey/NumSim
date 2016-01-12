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

#include "communicator.hpp"
#include "geometry.hpp"
#include "grid.hpp"
#include <mpi.h>
#include "iterator.hpp"
#include <cmath>

/**
 * Constuctor
 */
Communicator::Communicator(int* argc, char*** argv) {
	// Initialize MPI
	MPI_Init(argc,argv);

	MPI_Comm_size(MPI_COMM_WORLD,&_size);
	MPI_Comm_rank(MPI_COMM_WORLD,&_rank);

	printf("Hello from %d\n",_rank);
	printf("Numprocs is %d\n",_size);

	// Get max possible screen number
	index_t maxDivisor = (index_t)sqrt(_size);
	for ( index_t i = maxDivisor; i > 0; i--)
	{
		if (getSize() % i == 0)
		{
			maxDivisor = i;
			break;
		}
	}

	// Set rank rows and columns
	_tdim[0] = getSize() / maxDivisor;
	_tdim[1] = maxDivisor;

	//Sebbi
	_tidx[1] = (int)(getRank() / _tdim[0]);
	_tidx[0] = (int)(getRank() - _tidx[1] * _tdim[0]);

	//_tidx[0] = getRank() % maxDivisor;
	//_tidx[1] = (int)((getRank() - _tidx[0]) / maxDivisor);


	printf("Dims: %i %i \n", _tdim[0], _tdim[1]);
	printf("Rank: %i Idx : %i %i \n", getRank(), _tidx[0], _tidx[1]);

}

/**
 * Destructor
 */
Communicator::~Communicator() {
	MPI_Finalize();
}

/**
 * Returns the position of the current process with respect to the
 * fields lower left corner
 */
const multi_index_t& Communicator::ThreadIdx() const {
	return _tidx;
}

/**
 * Returns the way the domain is partitioned among all processes
 */
const multi_index_t& Communicator::ThreadDim() const {
	return _tdim;
}

/**
 * Returns whether this process is a red or a black field
 */
const bool& Communicator::EvenOdd() const {
	return _evenodd;
}

void Communicator::SetEvenOdd(bool evenodd) {
	_evenodd = evenodd;
}


/**
 * Gets the sum of all values and distributes the result among all
 * processes
 *
 * @param val  the data over which the sum is to be calculated
 */
real_t Communicator::gatherSum(const real_t& val) const {
	real_t valid = val, getResult;

	MPI_Allreduce(&valid, &getResult, 1, MPI_REAL_TYPE, MPI_SUM, MPI_COMM_WORLD);

	return getResult;
}

/**
 * Finds the minimum of the values and distributes the result among
 * all processes
 *
 * @param val  the data over which to find the minimum
 */
real_t Communicator::gatherMin(const real_t& val) const {
	real_t valid = val, getResult;

	MPI_Allreduce(&valid, &getResult, 1, MPI_REAL_TYPE, MPI_MIN, MPI_COMM_WORLD);

	return getResult;
}

/**
 * Finds the maximum of the values and distributes the result among
 * all processes
 *
 * @param val  the data over which to find the maximum
 */
real_t Communicator::gatherMax(const real_t& val) const {
	real_t valid = val, getResult;

	MPI_Allreduce(&valid, &getResult, 1, MPI_REAL_TYPE, MPI_MAX, MPI_COMM_WORLD);

	return getResult;
}

/**
 * Synchronizes ghost layer
 *
 * @param grid  the values to sync
 */
void Communicator::copyBoundary(Grid* grid) const {
	// Left, right
	if ( !(isLeft() && isRight() )) {
		copyLeftBoundary(grid);
		copyRightBoundary(grid);
	}

	// Top, bottom
	if ( ! (isTop() && isBottom()) ) {
		copyTopBoundary(grid);
		copyBottomBoundary(grid);
	}
}

/**
 * Decide whether our left boundary is a domain boundary
 */
bool Communicator::isLeft() const {
	//printf("L_curRak: %i \n", _rank);
	if ( ThreadIdx()[0] == 0 ) {
		//printf("yep \n");
		return true;
	}

	return false;
}

/**
 * Decide whether our right boundary is a domain boundary
 */
bool Communicator::isRight() const {
	//printf("R_curRak: %i %i %i \n", _rank, ThreadIdx()[0], ThreadDim()[0]);
	if ( (ThreadIdx()[0] + 1) == ThreadDim()[0] ) {
		//printf("yep \n");
		return true;
	}

	return false;
}

/**
 * Decide whether our top boundary is a domain boundary
 */
bool Communicator::isTop() const {
	//printf("T_curRak: %i \n", _rank);
	if ( (ThreadIdx()[1] + 1 ) == ThreadDim()[1] ) {
		//printf("yep \n");
		return true;
	}

	return false;
}

/**
 * Decide whether our bottom boundary is a domain boundary
 */
bool Communicator::isBottom() const {
	//printf("B_curRak: %i \n", _rank);
	if ( ThreadIdx()[1] == 0 ) {
		//printf("yep \n");
		return true;
	}

	return false;
}

/**
 * Get MPI rank of current process
 */
const int& Communicator::getRank() const {
	return _rank;
}

/**
 * Get number of MPI processes
 */
const int& Communicator::getSize() const {
	return _size;
}

/**
 * Function to search for the neighbour ranks
 *
 * @param side  get neighbor rank of rank side
 */
int Communicator::getNeighbour(int side) const {
	multi_index_t idx = _tidx;
	if      ( side == 0 ) {idx[1] -= 1;}
	else if ( side == 1 ) {idx[0] -= 1;}
	else if ( side == 2 ) {idx[1] += 1;}
	else if ( side == 3 ) {idx[0] += 1;}

	// rank = dim_x * idx_y + idx_x
	int rank = _tdim[0]*idx[1] + idx[0];
	//printf("curRak: %i dir: %i, rank: %i \n", _rank, side, rank);
	return rank;
}

/** Function to sync ghost layer on left boundary:
 *  send values of own left boundary to left neighbor and
 *  and receive values from his right boundary
 *
 * @param grid  values whose boundary shall be synced
 */
bool Communicator::copyLeftBoundary(Grid* grid) const {
	// Initialize
	real_t boundary[grid->getGeometry()->Size()[1]];
	real_t ghostLayer[grid->getGeometry()->Size()[1]];
	BoundaryIterator it = BoundaryIterator(grid->getGeometry());
	it.SetBoundary(it.boundaryLeft);
	int i = 0;
	// Cycle through boundary
	while(it.Valid()) {
		boundary[i] = grid->Cell(it.Right());
		it.Next();
		i++;
	}
	// If is right domain, update ghostlayer of neighbor
	if(isRight()) {
		//printf("left: send only \n");
		MPI_Send(&boundary,grid->getGeometry()->Size()[1], MPI_REAL_TYPE,
					getNeighbour((int)it.boundaryLeft), 0,MPI_COMM_WORLD);
	}
	// If is left domain, update own ghostlayer by neighbor
	else if(isLeft()){
		//printf("rank: %i \n", getRank());
		//printf("left: rec only \n");
		MPI_Recv(&ghostLayer, grid->getGeometry()->Size()[1], MPI_REAL_TYPE,
					getNeighbour((int)it.boundaryRight), 0,MPI_COMM_WORLD, MPI_STATUS_IGNORE);

	}
	// If is no domain, update left boundary of neighbor and update own right ghostlayer
	else{
		MPI_Sendrecv(&boundary,grid->getGeometry()->Size()[1], MPI_REAL_TYPE,
					getNeighbour((int)it.boundaryLeft), 0,
	                &ghostLayer, grid->getGeometry()->Size()[1], MPI_REAL_TYPE,
					getNeighbour((int)it.boundaryRight), 0,
					MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}
	// Set zeros for empty ghostlayer on first run
	if(!isRight()){
		it.SetBoundary(it.boundaryRight);
		i = 0;
		while(it.Valid()){
				grid->Cell(it) = ghostLayer[i];
				it.Next();
				i++;
		}
	}

	return true;
}

/**
 * Function to sync ghost layer on right boundary
 *  Details analog to left boundary
 *
 * @param grid  values whose boundary shall be synced
 */
bool Communicator::copyRightBoundary(Grid* grid) const {
	// Initialize
	real_t boundary[grid->getGeometry()->Size()[1]];
	real_t ghostLayer[grid->getGeometry()->Size()[1]];
	BoundaryIterator it = BoundaryIterator(grid->getGeometry());
	it.SetBoundary(it.boundaryRight);
	int i = 0;
	// Cycle through boundary
	while(it.Valid()){
		boundary[i] = grid->Cell(it.Left());
		it.Next();
		i++;
	}
	// If is left domain, update ghostlayer of neighbor
	if (isLeft()){
		MPI_Send(&boundary,grid->getGeometry()->Size()[1], MPI_REAL_TYPE,
					getNeighbour((int)it.boundaryRight), 0,MPI_COMM_WORLD);
	}
	// If is right domain, update own ghostlayer by neighbor
	else if(isRight()){
		int targetRank = getNeighbour((int)it.boundaryLeft);
		MPI_Recv(&ghostLayer, grid->getGeometry()->Size()[1], MPI_REAL_TYPE,
					targetRank , 0,MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}
	// If is no domain, update right boundary of neighbor and update own left ghostlayer
	else{
		MPI_Sendrecv(&boundary,grid->getGeometry()->Size()[1], MPI_REAL_TYPE,
					getNeighbour((int)it.boundaryRight), 0,
	                &ghostLayer, grid->getGeometry()->Size()[1], MPI_REAL_TYPE,
					getNeighbour((int)it.boundaryLeft), 0,
					MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}
	// Set zeros for empty ghostlayer on first run
	if(!isLeft()){
		it.SetBoundary(it.boundaryLeft);
		i = 0;
		while(it.Valid()){
				grid->Cell(it) = ghostLayer[i];
				it.Next();
				i++;
		}
	}

	return true;
}

/**
 * Function to sync ghost layer on top boundary
 *  Details analog to left boundary
 *
 * @param grid  values whose boundary shall be synced
 */
bool Communicator::copyTopBoundary(Grid* grid) const {
	// Initialize
	real_t boundary[grid->getGeometry()->Size()[0]];
	real_t ghostLayer[grid->getGeometry()->Size()[0]];
	BoundaryIterator it = BoundaryIterator(grid->getGeometry());
	it.SetBoundary(it.boundaryTop);
	int i = 0;
	// Cycle through boundary
	while(it.Valid()){
		boundary[i] = grid->Cell(it.Down());
		it.Next();
		i++;
	}
	// If is bottom domain, update ghostlayer of neighbor
	if (isBottom()){
		MPI_Send(&boundary,grid->getGeometry()->Size()[0], MPI_REAL_TYPE,
					getNeighbour((int)it.boundaryTop), 0,MPI_COMM_WORLD);
	}
	// If is top domain, update own ghostlayer by neighbor
	else if(isTop()){
		MPI_Recv(&ghostLayer, grid->getGeometry()->Size()[0], MPI_REAL_TYPE,
					getNeighbour((int)it.boundaryBottom), 0,MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}
	// If is no domain, update top boundary of neighbor and update own bottom ghostlayer
	else{
		MPI_Sendrecv(&boundary,grid->getGeometry()->Size()[0], MPI_REAL_TYPE,
					getNeighbour((int)it.boundaryTop), 0,
	                &ghostLayer, grid->getGeometry()->Size()[0], MPI_REAL_TYPE,
					getNeighbour((int)it.boundaryBottom), 0,
					MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}
	// Set zeros for empty ghostlayer on first run
	if(!isBottom()){
		it.SetBoundary(it.boundaryBottom);
		i = 0;
		while(it.Valid()){
			 	grid->Cell(it) = ghostLayer[i];
				it.Next();
				i++;
		}
	}

	return true;
}

/**
 * Function to sync ghost layer on bottom boundary
 * Details analog to left boundary
 *
 * @param grid  values whose boundary shall be synced
 */
bool Communicator::copyBottomBoundary(Grid* grid) const {
	// Initialize
	real_t boundary[grid->getGeometry()->Size()[0]];
	real_t ghostLayer[grid->getGeometry()->Size()[0]];
	BoundaryIterator it = BoundaryIterator(grid->getGeometry());
	it.SetBoundary(it.boundaryBottom);
	int i = 0;
	while(it.Valid()){
		boundary[i] = grid->Cell(it.Top());
		it.Next();
		i++;
	}
	// If is top domain, update ghostlayer of neighbor
	if (isTop()){
		MPI_Send(&boundary,grid->getGeometry()->Size()[0], MPI_REAL_TYPE,
					getNeighbour((int)it.boundaryBottom), 0,MPI_COMM_WORLD);
	}
	// If is bottom domain, update own ghostlayer by neighbor
	else if(isBottom()){
		MPI_Recv(&ghostLayer, grid->getGeometry()->Size()[0], MPI_REAL_TYPE,
					getNeighbour((int)it.boundaryTop), 0,MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}
	// If is no domain, update bottom boundary of neighbor and update own top ghostlayer
	else{
		MPI_Sendrecv(&boundary,grid->getGeometry()->Size()[0], MPI_REAL_TYPE,
					getNeighbour((int)it.boundaryBottom), 0,
	                &ghostLayer, grid->getGeometry()->Size()[0], MPI_REAL_TYPE,
					getNeighbour((int)it.boundaryTop), 0,
					MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}
	// Set zeros for empty ghostlayer on first run
	if(!isTop()){
		it.SetBoundary(it.boundaryTop);
		i = 0;
		while(it.Valid()){
			 	grid->Cell(it) = ghostLayer[i];
				it.Next();
				i++;
		}
	}

	return true;
}

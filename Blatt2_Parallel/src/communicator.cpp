#include "communicator.hpp"
#include "geometry.hpp"
#include "grid.hpp"
#include <mpich/mpi.h>
#include "iterator.hpp"
#include <cmath>

Communicator::Communicator(int* argc, char*** argv) {

	MPI_Init(argc,argv);

	MPI_Comm_size(MPI_COMM_WORLD,&_size);
	MPI_Comm_rank(MPI_COMM_WORLD,&_rank);

	printf("Hello from %d\n",_rank);
	printf("Numprocs is %d\n",_size);


	index_t maxDivisor = (index_t)sqrt(_size);
	for ( index_t i = maxDivisor; i > 0; i--)
	{
		if (getSize() % i == 0)
		{
			maxDivisor = i;
			break;
		}
	}

	_tdim[0] = getSize() / maxDivisor;
	_tdim[1] = maxDivisor;

	_tidx[0] = (int)(getRank() / maxDivisor);
	_tidx[1] = getRank() % maxDivisor;

	printf("Dims: %i %i \n", _tdim[0], _tdim[1]);
	printf("Idx : %i %i \n", _tidx[0], _tidx[1]);

}

Communicator::~Communicator() {
	MPI_Finalize();
}

const multi_index_t& Communicator::ThreadIdx() const {
	return _tidx;
}

const multi_index_t& Communicator::ThreadDim() const {
	return _tdim;
}

const bool& Communicator::EvenOdd() const {
	return _evenodd;
}

real_t Communicator::gatherSum(const real_t& val) const {
	double valid = val, getResult;

	MPI_Allreduce(&valid, &getResult, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

	return getResult;
}

real_t Communicator::gatherMin(const real_t& val) const {
	double valid = val, getResult;

	MPI_Allreduce(&valid, &getResult, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

	return getResult;
}

real_t Communicator::gatherMax(const real_t& val) const {
	double valid = val, getResult;

	MPI_Allreduce(&valid, &getResult, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

	return getResult;
}

void Communicator::copyBoundary(Grid* grid) const {
	copyLeftBoundary(grid);
	copyRightBoundary(grid);
	copyTopBoundary(grid);
	copyBottomBoundary(grid);
}

const bool Communicator::isLeft() const {
}

const bool Communicator::isRight() const {
}

const bool Communicator::isTop() const {
}

const bool Communicator::isBottom() const {
}

const int& Communicator::getRank() const {
	return _rank;
}

const int& Communicator::getSize() const {
	return _size;
}

int Communicator::getNeighbour(int side) const {
	multi_index_t idx = _tidx;
	if ( side == 0 )
	{
		idx[1] -= 1;
	}
	else if ( side == 1 )
	{
		idx[0] -= 1;
	}
	else if ( side == 2 )
	{
		idx[1] += 1;
	}
	else if ( side == 3 )
	{
		idx[0] += 1;
	}

	// rank = dim_x * idx_y + idx_x
	return _tdim[1]*idx[0] + idx[1];
}

bool Communicator::copyLeftBoundary(Grid* grid) const {
	real_t boundary[grid->getGeometry()->Size()[1]];
	real_t ghostLayer[grid->getGeometry()->Size()[1]];
	BoundaryIterator it = BoundaryIterator(grid->getGeometry());
	it.SetBoundary(it.boundaryLeft);
	int i = 0;
	while(it.Valid()){
		boundary[i] = grid->Cell(it.Right());
		it.Next();
		i++;
	}
	if (isRight()){
		MPI_Send(&boundary,grid->getGeometry()->Size()[1], MPI_Datatype MPI_REAL_TYPE,
					getNeighbour((int)it.boundaryLeft), 0,MPI_COMM_WORLD);
	}
	else if(isLeft()){
		MPI_Recv(&ghostLayer, grid->getGeometry()->Size()[1], MPI_Datatype MPI_REAL_TYPE,
					getNeighbour((int)it.boundaryRight), 0,MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}
	else{
		MPI_Sendrecv(&boundary,grid->getGeometry()->Size()[1], MPI_Datatype MPI_REAL_TYPE,
					getNeighbour((int)it.boundaryLeft), 0,
	                &ghostLayer, grid->getGeometry()->Size()[1], MPI_Datatype MPI_REAL_TYPE,
					getNeighbour((int)it.boundaryRight), 0,
					MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}


	it.SetBoundary(it.boundaryRight);
	i = 0;
	while(it.Valid()){
		 	grid->Cell(it) = ghostLayer[i];
			it.Next();
			i++;
	}
}

bool Communicator::copyRightBoundary(Grid* grid) const {
	real_t boundary[grid->getGeometry()->Size()[1]];
	real_t ghostLayer[grid->getGeometry()->Size()[1]];
	BoundaryIterator it = BoundaryIterator(grid->getGeometry());
	it.SetBoundary(it.boundaryRight);
	int i = 0;
	while(it.Valid()){
		boundary[i] = grid->Cell(it.Left());
		it.Next();
		i++;
	}
	if (isLeft()){
		MPI_Send(&boundary,grid->getGeometry()->Size()[1], MPI_Datatype MPI_REAL_TYPE,
					getNeighbour((int)it.boundaryRight), 0,MPI_COMM_WORLD);
	}
	else if(isRight()){
		MPI_Recv(&ghostLayer, grid->getGeometry()->Size()[1], MPI_Datatype MPI_REAL_TYPE,
					getNeighbour((int)it.boundaryLeft), 0,MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}
	else{
		MPI_Sendrecv(&boundary,grid->getGeometry()->Size()[1], MPI_Datatype MPI_REAL_TYPE,
					getNeighbour((int)it.boundaryRight), 0,
	                &ghostLayer, grid->getGeometry()->Size()[1], MPI_Datatype MPI_REAL_TYPE,
					getNeighbour((int)it.boundaryLeft), 0,
					MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}


	it.SetBoundary(it.boundaryLeft);
	i = 0;
	while(it.Valid()){
		 	grid->Cell(it) = ghostLayer[i];
			it.Next();
			i++;
	}
}

bool Communicator::copyTopBoundary(Grid* grid) const {
	real_t boundary[grid->getGeometry()->Size()[1]];
	real_t ghostLayer[grid->getGeometry()->Size()[1]];
	BoundaryIterator it = BoundaryIterator(grid->getGeometry());
	it.SetBoundary(it.boundaryTop);
	int i = 0;
	while(it.Valid()){
		boundary[i] = grid->Cell(it.Down());
		it.Next();
		i++;
	}
	if (isBottom()){
		MPI_Send(&boundary,grid->getGeometry()->Size()[1], MPI_Datatype MPI_REAL_TYPE,
					getNeighbour((int)it.boundaryTop), 0,MPI_COMM_WORLD);
	}
	else if(isTop()){
		MPI_Recv(&ghostLayer, grid->getGeometry()->Size()[1], MPI_Datatype MPI_REAL_TYPE,
					getNeighbour((int)it.boundaryBottom), 0,MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}
	else{
		MPI_Sendrecv(&boundary,grid->getGeometry()->Size()[1], MPI_Datatype MPI_REAL_TYPE,
					getNeighbour((int)it.boundaryTop), 0,
	                &ghostLayer, grid->getGeometry()->Size()[1], MPI_Datatype MPI_REAL_TYPE,
					getNeighbour((int)it.boundaryBottom), 0,
					MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}


	it.SetBoundary(it.boundaryBottom);
	i = 0;
	while(it.Valid()){
		 	grid->Cell(it) = ghostLayer[i];
			it.Next();
			i++;
	}
}

bool Communicator::copyBottomBoundary(Grid* grid) const {
	real_t boundary[grid->getGeometry()->Size()[1]];
	real_t ghostLayer[grid->getGeometry()->Size()[1]];
	BoundaryIterator it = BoundaryIterator(grid->getGeometry());
	it.SetBoundary(it.boundaryBottom);
	int i = 0;
	while(it.Valid()){
		boundary[i] = grid->Cell(it.Top());
		it.Next();
		i++;
	}
	if (isTop()){
		MPI_Send(&boundary,grid->getGeometry()->Size()[1], MPI_Datatype MPI_REAL_TYPE,
					getNeighbour((int)it.boundaryBottom), 0,MPI_COMM_WORLD);
	}
	else if(isBottom()){
		MPI_Recv(&ghostLayer, grid->getGeometry()->Size()[1], MPI_Datatype MPI_REAL_TYPE,
					getNeighbour((int)it.boundaryTop), 0,MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}
	else{
		MPI_Sendrecv(&boundary,grid->getGeometry()->Size()[1], MPI_Datatype MPI_REAL_TYPE,
					getNeighbour((int)it.boundaryBottom), 0,
	                &ghostLayer, grid->getGeometry()->Size()[1], MPI_Datatype MPI_REAL_TYPE,
					getNeighbour((int)it.boundaryTop), 0,
					MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}


	it.SetBoundary(it.boundaryTop);
	i = 0;
	while(it.Valid()){
		 	grid->Cell(it) = ghostLayer[i];
			it.Next();
			i++;
	}
}

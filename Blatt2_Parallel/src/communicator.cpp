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
}

real_t Communicator::gatherMin(const real_t& val) const {
}

real_t Communicator::gatherMax(const real_t& val) const {
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

int Communicator::getNeighbour(int side) {
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
	BoundaryIterator it = BoundaryIterator(grid->getGeometry());
	it.SetBoundary(it.boundaryLeft);
}

bool Communicator::copyRightBoundary(Grid* grid) const {
}

bool Communicator::copyTopBoundary(Grid* grid) const {
}

bool Communicator::copyBottomBoundary(Grid* grid) const {
}

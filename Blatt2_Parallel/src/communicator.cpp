#include "communicator.hpp"
#include "geometry.hpp"
#include "grid.hpp"
#include <mpich/mpi.h>
#include "iterator.hpp"

Communicator::Communicator(int* argc, char*** argv) {

	MPI_Init(argc,argv);
	MPI_Comm_size(MPI_COMM_WORLD,&_size);
	MPI_Comm_rank(MPI_COMM_WORLD,&_rank);

	printf("Hello from %d\n",_rank);
	printf("Numprocs is %d\n",_size);

	_tidx[0] = (**argv)[0];
	_tidx[1] = (**argv)[1];

	_tdim[0] = (**argv)[2];
	_tdim[1] = (**argv)[3];

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
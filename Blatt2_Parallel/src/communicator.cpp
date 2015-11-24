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

bool Communicator::copyLeftBoundary(Grid* grid) const {
	real_t leftBoundary[grid->getGeometry()->Size()[1]];
	real_t ghostLayer[grid->getGeometry()->Size()[1]];

	//leftBoundary = new real_t [grid->_geom->Size()[1]];
	BoundaryIterator it = BoundaryIterator(grid->getGeometry());
	it.SetBoundary(it.boundaryLeft);
	int i = 0;
	while(it.Valid()){
		leftBoundary[i] = grid->Cell(it);
		it.Next();
		i++;
	}

	int b = 2;
	MPI_Sendrecv(&leftBoundary,grid->getGeometry()->Size()[1], MPI_Datatype MPI_REAL_TYPE,
					getNeighbour((int)it.boundaryLeft), 0,
	                &ghostLayer, grid->getGeometry()->Size()[1], MPI_Datatype MPI_REAL_TYPE,
					getNeighbour((int)it.boundaryRight), 0,
					MPI_COMM_WORLD, MPI_STATUS_IGNORE);

	it.SetBoundary(it.boundaryRight);
	i = 0;
	while(it.Valid()){
		 	grid->Cell(it.Left()) = leftBoundary[i];
			it.Next();
			i++;
	}
	int a = getNeighbour((int)it.boundaryRight);

}

bool Communicator::copyRightBoundary(Grid* grid) const {
}

bool Communicator::copyTopBoundary(Grid* grid) const {
}

bool Communicator::copyBottomBoundary(Grid* grid) const {
}

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

	//Sebbi
	_tidx[1] = (int)(getRank() / _tdim[0]);
	_tidx[0] = (int)(getRank() - _tidx[1] * _tdim[0]);

	//_tidx[0] = getRank() % maxDivisor;
	//_tidx[1] = (int)((getRank() - _tidx[0]) / maxDivisor);


	printf("Dims: %i %i \n", _tdim[0], _tdim[1]);
	printf("Rank: %i Idx : %i %i \n", getRank(), _tidx[0], _tidx[1]);

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

void Communicator::SetEvenOdd(bool evenodd) {
	_evenodd = evenodd;
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

bool Communicator::isLeft() const {
	//printf("L_curRak: %i \n", _rank);
	if ( ThreadIdx()[0] == 0 )
	{
		//printf("yep \n");
		return true;
	}

	return false;
}

bool Communicator::isRight() const {
	//printf("R_curRak: %i %i %i \n", _rank, ThreadIdx()[0], ThreadDim()[0]);
	if ( (ThreadIdx()[0] + 1) == ThreadDim()[0] )
	{
		//printf("yep \n");
		return true;
	}

	return false;
}

bool Communicator::isTop() const {
	//printf("T_curRak: %i \n", _rank);
	if ( (ThreadIdx()[1] + 1 ) == ThreadDim()[1] )
	{
		//printf("yep \n");
		return true;
	}



	return false;
}

bool Communicator::isBottom() const {
	//printf("B_curRak: %i \n", _rank);
	if ( ThreadIdx()[1] == 0 )
	{
		//printf("yep \n");
		return true;
	}

	return false;
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
	int rank = _tdim[0]*idx[1] + idx[0];
	//printf("curRak: %i dir: %i, rank: %i \n", _rank, side, rank);
	return rank;
}

bool Communicator::copyLeftBoundary(Grid* grid) const {
	if ( isRight() && isLeft())
		return true;

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
	if (isRight() && !isLeft()){
		//printf("left: send only \n");
		MPI_Send(&boundary,grid->getGeometry()->Size()[1], MPI_Datatype MPI_REAL_TYPE,
					getNeighbour((int)it.boundaryLeft), 0,MPI_COMM_WORLD);

	}
	else if(isLeft() && !isRight()){
		//printf("rank: %i \n", getRank());
		//printf("left: rec only \n");
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
	return true;
}

bool Communicator::copyRightBoundary(Grid* grid) const {
	if ( isRight() && isLeft())
			return true;

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
	if (isLeft() && !isRight()){
		MPI_Send(&boundary,grid->getGeometry()->Size()[1], MPI_Datatype MPI_REAL_TYPE,
					getNeighbour((int)it.boundaryRight), 0,MPI_COMM_WORLD);
	}
	else if(isRight() && !isLeft()){
		int targetRank = getNeighbour((int)it.boundaryLeft);
		MPI_Recv(&ghostLayer, grid->getGeometry()->Size()[1], MPI_Datatype MPI_REAL_TYPE,
					targetRank , 0,MPI_COMM_WORLD, MPI_STATUS_IGNORE);
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
	return true;
}

bool Communicator::copyTopBoundary(Grid* grid) const {
	if ( isTop() && isBottom())
			return true;
	real_t boundary[grid->getGeometry()->Size()[0]];
	real_t ghostLayer[grid->getGeometry()->Size()[0]];
	BoundaryIterator it = BoundaryIterator(grid->getGeometry());
	it.SetBoundary(it.boundaryTop);
	int i = 0;
	while(it.Valid()){
		boundary[i] = grid->Cell(it.Down());
		it.Next();
		i++;
	}
	if (isBottom() && !isTop()){
		MPI_Send(&boundary,grid->getGeometry()->Size()[0], MPI_Datatype MPI_REAL_TYPE,
					getNeighbour((int)it.boundaryTop), 0,MPI_COMM_WORLD);
	}
	else if(isTop() && !isBottom()){
		MPI_Recv(&ghostLayer, grid->getGeometry()->Size()[0], MPI_Datatype MPI_REAL_TYPE,
					getNeighbour((int)it.boundaryBottom), 0,MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}
	else{
		MPI_Sendrecv(&boundary,grid->getGeometry()->Size()[0], MPI_Datatype MPI_REAL_TYPE,
					getNeighbour((int)it.boundaryTop), 0,
	                &ghostLayer, grid->getGeometry()->Size()[0], MPI_Datatype MPI_REAL_TYPE,
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
	return true;
}

bool Communicator::copyBottomBoundary(Grid* grid) const {
	if ( isTop() && isBottom())
		return true;

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
	if (isTop() && !isBottom()){
		MPI_Send(&boundary,grid->getGeometry()->Size()[0], MPI_Datatype MPI_REAL_TYPE,
					getNeighbour((int)it.boundaryBottom), 0,MPI_COMM_WORLD);
	}
	else if(isBottom() && !isTop()){
		MPI_Recv(&ghostLayer, grid->getGeometry()->Size()[0], MPI_Datatype MPI_REAL_TYPE,
					getNeighbour((int)it.boundaryTop), 0,MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}
	else{
		MPI_Sendrecv(&boundary,grid->getGeometry()->Size()[0], MPI_Datatype MPI_REAL_TYPE,
					getNeighbour((int)it.boundaryBottom), 0,
	                &ghostLayer, grid->getGeometry()->Size()[0], MPI_Datatype MPI_REAL_TYPE,
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
	return true;
}

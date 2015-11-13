#include "communicator.hpp"
#include "geometry.hpp"
#include "grid.hpp"

Communicator::Communicator(int* argc, char*** argv) {
}

Communicator::~Communicator() {
}

const multi_index_t& Communicator::ThreadIdx() const {
}

const multi_index_t& Communicator::ThreadDim() const {
}

const bool& Communicator::EvenOdd() const {
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
}

const int& Communicator::getSize() const {
}

bool Communicator::copyLeftBoundary(Grid* grid) const {
}

bool Communicator::copyRightBoundary(Grid* grid) const {
}

bool Communicator::copyTopBoundary(Grid* grid) const {
}

bool Communicator::copyBottomBoundary(Grid* grid) const {
}

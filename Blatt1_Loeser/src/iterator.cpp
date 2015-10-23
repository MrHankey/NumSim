#include "iterator.hpp"

Iterator::Iterator(const Geometry* geom) {
}

Iterator::Iterator(const Geometry* geom, const index_t& value) {
}

const index_t& Iterator::Value() const {
}

Iterator::operator const index_t&() const {
}

multi_index_t Iterator::Pos() const {
}

void Iterator::First() {
}

void Iterator::Next() {
}

bool Iterator::Valid() const {
}

Iterator Iterator::Left() const {
}

Iterator Iterator::Right() const {
}

Iterator Iterator::Top() const {
}

Iterator Iterator::Down() const {
}

InteriorIterator::InteriorIterator(const Geometry* geom) : Iterator(geom) {

}

void InteriorIterator::First() {
}

void InteriorIterator::Next() {
}

BoundaryIterator::BoundaryIterator(const Geometry* geom) : Iterator(geom) {
}

void BoundaryIterator::SetBoundary(const index_t& boundary) {
}

void BoundaryIterator::First() {
}

void BoundaryIterator::Next() {
}

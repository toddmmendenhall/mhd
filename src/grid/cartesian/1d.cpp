#include <1d.hpp>
#include <grid.hpp>
#include <profile.hpp>

#include <memory>

namespace MHD {

Cartesian1DGrid::Cartesian1DGrid(Profile const& profile) : Grid(profile) {}

Cartesian1DGrid::~Cartesian1DGrid() = default;

void Cartesian1DGrid::SomeMethod() {}

} // namespace MHD

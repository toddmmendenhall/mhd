#include <3d.hpp>
#include <grid.hpp>
#include <profile.hpp>

#include <memory>

namespace MHD {

Cartesian3DGrid::Cartesian3DGrid(Profile const& profile) : Grid(profile) {}

Cartesian3DGrid::~Cartesian3DGrid() = default;

void Cartesian3DGrid::SomeMethod() {}

} // namespace MHD

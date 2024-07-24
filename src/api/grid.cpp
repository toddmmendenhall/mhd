#include "grid.hpp"
#include "grid_factory.hpp"
#include "grid_impl.hpp"
#include "profile.hpp"

#include <memory>

namespace MHD {

Grid::Grid(std::unique_ptr<Profile> const& profile) {
    GridFactory gridFactory;
    m_gridImpl = gridFactory.CreateGrid(profile);
}

Grid::~Grid() {}

} // namespace MHD

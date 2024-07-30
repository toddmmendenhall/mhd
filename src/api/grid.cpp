#include "grid.hpp"
#include "grid_factory.hpp"
#include "grid_impl.hpp"
#include "profile.hpp"

#include <memory>

namespace MHD {

Grid::Grid(std::unique_ptr<Profile> const& profile) {
    auto gridFactory = new GridFactory();
    m_gridImpl = gridFactory->CreateGrid(profile);
    delete gridFactory;
}

Grid::~Grid() {}

std::unique_ptr<GridImpl> const& Grid::GetGridImpl() const {
    return m_gridImpl;
}

} // namespace MHD

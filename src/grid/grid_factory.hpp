#pragma once

#include "grid_impl.hpp"
#include "profile.hpp"

#include <memory>

namespace MHD {

class GridFactory {
public:
    std::unique_ptr<GridImpl> CreateGrid(std::unique_ptr<Profile> const& profile) const;
};

} // namespace MHD

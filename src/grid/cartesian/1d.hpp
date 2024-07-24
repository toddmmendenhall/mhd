#pragma once

#include "grid_impl.hpp"
#include "profile.hpp"

#include <memory>

namespace MHD {

class Cartesian1DGrid : public GridImpl {
public:
    Cartesian1DGrid(std::unique_ptr<Profile> const& profile);
};

} // namespace MHD

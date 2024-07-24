#pragma once

#include "grid_impl.hpp"
#include "profile.hpp"

#include <memory>

namespace MHD {

class Cartesian3DGrid : public GridImpl {
public:
    Cartesian3DGrid(std::unique_ptr<Profile> const& profile);
};

} // namespace MHD

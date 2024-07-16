#pragma once

#include "grid.hpp"
#include "profile.hpp"

namespace MHD {

class GridFactory {
public:
    Grid* CreateGrid(Profile* profile) const;
};

} // namespace MHD

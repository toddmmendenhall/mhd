#pragma once

#include <grid.hpp>

#include <memory>
#include <vector>

namespace MHD {

class Profile;

class Cartesian1DGrid : public IGrid {
public:
    Cartesian1DGrid(Profile const& profile);
};

} // namespace MHD

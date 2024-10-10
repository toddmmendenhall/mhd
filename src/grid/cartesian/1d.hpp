#pragma once

#include <grid.hpp>

#include <memory>

namespace MHD {

class Profile;

class Cartesian1DGrid : public Grid {
public:
    Cartesian1DGrid(Profile const& profile);
    ~Cartesian1DGrid();

    void SomeMethod() override;
};

} // namespace MHD

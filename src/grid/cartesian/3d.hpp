#pragma once

#include <grid.hpp>

#include <memory>

namespace MHD {

class Profile;

class Cartesian3DGrid : public IGrid {
public:
    Cartesian3DGrid(Profile const& profile);
    ~Cartesian3DGrid();

    void SomeMethod() override;
};

} // namespace MHD

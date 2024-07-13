#pragma once

#include "grid.hpp"
#include "profile.hpp"

namespace MHD {

class Calc {
public:
    Calc();
    ~Calc();

    Profile* GetProfile() const;
    Grid* GetGrid() const;

private:
    Profile* m_profile;
    Grid* m_grid;
};

} // namespace MHD

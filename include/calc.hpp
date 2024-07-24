#pragma once

#include "grid.hpp"
#include "profile.hpp"

#include <memory>

namespace MHD {

class Calc {
public:
    Calc();

    void Run();

    std::unique_ptr<Profile> const& GetProfile() const;
    std::unique_ptr<Grid> const& GetGrid() const;

private:
    void SetupCalc();

    std::unique_ptr<Profile> m_profile;
    std::unique_ptr<Grid> m_grid;
};

} // namespace MHD

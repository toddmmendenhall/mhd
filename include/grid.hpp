#pragma once

#include "profile.hpp"

#include <memory>

namespace MHD {

class GridImpl;

class Grid {
public:
    Grid(std::unique_ptr<Profile> const& profile);
    ~Grid();

    std::unique_ptr<GridImpl> const& GetGridImpl() const;

private:
    std::unique_ptr<GridImpl> m_gridImpl;
};

} // namespace MHD

#pragma once

#include "grid.hpp"
#include "profile.hpp"

#include <memory>
#include <vector>

namespace MHD {

class CartesianGrid : public Grid {
public:
    CartesianGrid(std::unique_ptr<GridProfile> const& gridProfile);
    ~CartesianGrid();

private:
    /// @brief List of nodes in domain
    std::vector<std::vector<double>> m_nodes;
};

} // namespace MHD

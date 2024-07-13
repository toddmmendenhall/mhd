#pragma once

#include "grid.hpp"
#include "profile_options.hpp"

#include <array>
#include <vector>

namespace MHD {

class Cartesian1DGrid : public Grid {
public:
    Cartesian1DGrid(std::vector<double> const& bounds, std::vector<double> const& spacing, std::vector<BoundaryCondition> const& boundaryConditions);
};

class Cartesian2DGrid : public Grid {
public:
    Cartesian2DGrid(std::vector<double> const& bounds, std::vector<double> const& spacing, std::vector<BoundaryCondition> const& boundaryConditions);
private:
    double m_xMin;
    double m_xMax;
    double m_yMin;
    double m_yMax;
    std::vector<std::array<double, 2>> m_cells;
    std::vector<std::array<double, 2>> m_boundaryCells;
};

class Cartesian3DGrid : public Grid {
public:
    Cartesian3DGrid(std::vector<double> const& bounds, std::vector<double> const& spacing, std::vector<BoundaryCondition> const& boundaryConditions);
};

} // namespace MHD

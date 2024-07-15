#pragma once

#include "cell.hpp"
#include "grid.hpp"
#include "profile.hpp"
#include "profile_options.hpp"

#include <array>
#include <vector>

namespace MHD {

class Cartesian1DGrid : public Grid {
public:
    Cartesian1DGrid(Profile* profile);
};

class Cartesian2DGrid : public Grid {
public:
    Cartesian2DGrid(Profile* profile);

private:
    double m_xMin;
    double m_xMax;
    double m_yMin;
    double m_yMax;
    std::vector<Cell> m_interiorCells;
    std::vector<std::array<double, 2>> m_boundaryCells;
};

class Cartesian3DGrid : public Grid {
public:
    Cartesian3DGrid(Profile* profile);
};

} // namespace MHD

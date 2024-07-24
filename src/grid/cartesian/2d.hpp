#pragma once

#include "cell.hpp"
#include "grid_impl.hpp"
#include "profile.hpp"
#include "profile_options.hpp"

#include <array>
#include <memory>
#include <vector>

namespace MHD {

class Cartesian2DGrid : public GridImpl {
public:
    Cartesian2DGrid(std::unique_ptr<Profile> const& profile);

private:
    double m_xMin;
    double m_xMax;
    double m_yMin;
    double m_yMax;
    std::vector<Cell> m_interiorCells;
    std::vector<std::array<double, 2>> m_boundaryCells;
};

} // namespace MHD

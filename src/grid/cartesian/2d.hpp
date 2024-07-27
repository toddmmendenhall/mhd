#pragma once

#include "grid_impl.hpp"
#include "profile.hpp"
#include "point.hpp"
#include "vector.hpp"

#include <memory>
#include <vector>

namespace MHD {

class Cartesian2DGrid : public GridImpl {
public:
    Cartesian2DGrid(std::unique_ptr<Profile> const& profile);

private:
    std::vector<GEOM_UTILS::Point2D> m_cells;
    std::vector<bool> m_boundaryCellFlag;
    std::vector<double> m_density;
    std::vector<GEOM_UTILS::Vector2D> m_velocity;
    std::vector<double> m_pressure;
    std::vector<GEOM_UTILS::Vector2D> m_magneticField;
};

} // namespace MHD

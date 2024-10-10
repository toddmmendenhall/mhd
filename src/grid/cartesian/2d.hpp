#pragma once

#include <grid.hpp>
// #include <point.hpp>
// #include <vector.hpp>

#include <memory>
#include <vector>

namespace MHD {

class Profile;

class Cartesian2DGrid : public Grid {
public:
    Cartesian2DGrid(Profile const& profile);
    ~Cartesian2DGrid();

    void SomeMethod() override;

private:
    // std::vector<GEOM_UTILS::Point2D> m_cells;
    // std::vector<bool> m_boundaryCellFlag;
    // std::vector<double> m_density;
    // std::vector<GEOM_UTILS::Vector2D> m_velocity;
    // std::vector<double> m_pressure;
    // std::vector<GEOM_UTILS::Vector2D> m_magneticField;
};

} // namespace MHD

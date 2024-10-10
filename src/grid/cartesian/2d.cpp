#include <2d.hpp>
#include <grid.hpp>
// #include <constants.hpp>
// #include <point.hpp>

#include <cmath>
#include <memory>
#include <vector>

namespace MHD {

Cartesian2DGrid::Cartesian2DGrid(Profile const& profile) : Grid(profile) {
    // std::vector<double> const& bounds = profile->GetGridBounds();
    // std::vector<double> const& spacing = profile->GetGridSpacing();

    // double const xMin = bounds[0];
    // double const xMax = bounds[1];
    // double const yMin = bounds[2];
    // double const yMax = bounds[3];

    // // Calculates the closest integer
    // double const xNumCells = std::round((xMax - xMin) / spacing[0]);
    // double const yNumCells = std::round((yMax - yMin) / spacing[1]);

    // // Changes the user-defined spacings so that cells are equally distributed inside
    // // the domain
    // double const xSpacing = (xMax - xMin) / xNumCells;
    // double const ySpacing = (yMax - yMin) / yNumCells;

    // // The total number of cells also includes boundary cells that are not inside the domain
    // std::size_t const iMax = xNumCells + 2;
    // std::size_t const jMax = yNumCells + 2;
    // std::size_t const numCells = iMax * jMax;

    // m_cells.reserve(numCells);
    // m_boundaryCellFlag.assign(numCells, false);

    // // Set up cells and boundary cells
    // for (std::size_t i = 0; i < iMax; ++i) {
    //     for (std::size_t j = 0; j < jMax; ++j) {
    //         m_cells.push_back(GEOM_UTILS::Point2D(xMin + (i + 0.5) * xSpacing,
    //                                               yMin + (j + 0.5) * ySpacing));

    //         if (i == 0 || i == iMax - 1 || j == 0 || j == jMax - 1) {
    //             m_boundaryCellFlag[i * jMax + j] = true;
    //         }
    //     }
    // }

    // m_density.assign(numCells, ATMOSPHERIC_DENSITY_STP);
    // m_velocity.assign(numCells, GEOM_UTILS::Vector2D(0.0, 0.0));
    // m_pressure.assign(numCells, ATMOSPHERIC_PRESSURE_STP);
    // m_magneticField.assign(numCells, GEOM_UTILS::Vector2D(0.0, 0.0));
}

Cartesian2DGrid::~Cartesian2DGrid() = default;

void Cartesian2DGrid::SomeMethod() {}

} // namespace MHD

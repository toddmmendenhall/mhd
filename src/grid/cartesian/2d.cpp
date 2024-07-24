#include "2d.hpp"
#include "cell.hpp"

#include <cmath>
#include <memory>
#include <vector>

namespace MHD {



Cartesian2DGrid::Cartesian2DGrid(std::unique_ptr<Profile> const& profile) {
    std::vector<double> const& bounds = profile->GetGridBounds();
    std::vector<double> const& spacing = profile->GetGridSpacing();
    std::vector<BoundaryCondition> const& boundaryConditions = profile->GetGridBoundaryConditions();

    m_xMin = bounds[0];
    m_xMax = bounds[1];
    m_yMin = bounds[2];
    m_yMax = bounds[3];

    std::size_t const numCellsX = std::round((m_xMax - m_xMin) / spacing[0]);
    std::size_t const numCellsY = std::round((m_yMax - m_yMin) / spacing[1]);

    double const spacingX = (m_xMax - m_xMin) / numCellsX;
    double const spacingY = (m_yMax - m_yMin) / numCellsY;

    std::vector<std::vector<double>> cellPositions;
    for (std::size_t i = 0; i < numCellsX; ++i) {
        for (std::size_t j = 0; j < numCellsY; ++j) {
            double cellX = m_xMin + i * spacingX;
            double cellY = m_yMin + j * spacingY;

            cellPositions.push_back({cellX, cellY});

            if (j == 0) {
                m_boundaryCells.push_back({cellX, cellY - spacingY});
            }
            if (i == 0) {
                m_boundaryCells.push_back({cellX - spacingX, cellY});
            }
            if (i == numCellsX - 1) {
                m_boundaryCells.push_back({cellX + spacingX, cellY});
            }
            if (j == numCellsY - 1) {
                m_boundaryCells.push_back({cellX, cellY + spacingY});
            }
        }
    }
    //CellFactory* cellFactory = new CellFactory();
    //m_interiorCells = cellFactory->SetCells(profile, cellPositions);
    //delete cellFactory;
}



} // namespace MHD

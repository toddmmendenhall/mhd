#pragma once

#include "profile.hpp"
#include "profile_options.hpp"

#include <vector>

namespace MHD {

class Cell {};

class CellFactory {
public:
    std::vector<Cell> SetCells(Profile* profile, std::vector<std::vector<double>> const& positions) const {}
};

class FiniteDifferenceCell : public Cell {
public:
    FiniteDifferenceCell();

private:
    std::size_t m_index;
    std::vector<double> m_coordinates;
    std::vector<Cell*> m_nearestNeighbors;
    BoundaryCondition m_boundaryCondition;
    std::vector<double> m_state;
};

} // namespace MHD

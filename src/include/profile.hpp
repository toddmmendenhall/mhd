#pragma once

#include "profile_options.hpp"

#include <vector>

namespace MHD {

class Profile {
public:
    Profile();

    Dimension const GetGridDimension() const;
    Geometry const GetGridGeometry() const;
    std::vector<double> const& GetGridBounds() const;
    std::vector<double> const& GetGridSpacing() const;
    std::vector<BoundaryCondition> const& GetGridBoundaryConditions() const;

    void SetGridDimension(Dimension const dimension);
    void SetGridGeometry(Geometry const geometry);
    void SetGridBounds(std::vector<double> const& bounds);
    void SetGridSpacing(std::vector<double> const& spacing);
    void SetGridBoundaryConditions(std::vector<BoundaryCondition> const& boundaryConditions);

private:
    Dimension m_gridDimension;
    Geometry m_gridGeometry;
    std::vector<double> m_gridBounds;
    std::vector<double> m_gridSpacing;
    std::vector<BoundaryCondition> m_gridBoundaryConditions;
};

} // namespace MHD

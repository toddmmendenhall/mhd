#pragma once

#include <array>
#include <memory>

namespace MHD {

enum class Geometry {
    CARTESIAN = 0,
};

enum class BoundaryCondition {
    DIRICHLET = 0,
};

struct GridProfile {
    short m_dimension = 2;

    Geometry m_geometry = Geometry::CARTESIAN;

    std::array<BoundaryCondition, 6> m_boundaryConditions = {
        BoundaryCondition::DIRICHLET, // First coordinate lower boundary condition
        BoundaryCondition::DIRICHLET, // First coordinate upper boundary condition
        BoundaryCondition::DIRICHLET, // Second coordinate lower boundary condition
        BoundaryCondition::DIRICHLET, // Second coordinate upper boundary condition
        BoundaryCondition::DIRICHLET, // Third coordinate lower boundary condition
        BoundaryCondition::DIRICHLET  // Third coordinate upper boundary condition
    };

    std::array<double, 6> m_limits = {
        -5., // First coordinate minimum value
        5.,  // First coordinate maximum value
        -5., // Second coordinate minimum value
        5.,  // Second coordinate maximum value
        -5., // Third coordinate minimum value
        5.   // Third coordinate maximum value
    };
};

class Profile {
public:
    Profile();
    ~Profile();

    short const GetGridDimension() const;
    Geometry const GetGridGeometry() const;
    std::array<BoundaryCondition, 6> const& GetGridBoundaryConditions() const;
    std::array<double, 6> const& GetGridLimits() const;

    void SetGridDimension(short const dimension);
    void SetGridGeometry(Geometry const geometry);
    void SetGridBoundaryConditions(std::array<BoundaryCondition, 6> const& boundaryConditions);
    void SetGridLimits(std::array<double, 6> const& limits);

private:
    std::unique_ptr<GridProfile> m_gridProfile;
};

}

#pragma once

#include <vector>
#include <memory>

namespace MHD {

enum class Geometry {
    CARTESIAN = 0,
};

enum class BoundaryCondition {
    DIRICHLET = 0,
};

struct GridProfile {
    Geometry m_geometry = Geometry::CARTESIAN;

    unsigned short m_dimension = 2;

    std::vector<double> m_limits = {
        -5., // First coordinate minimum value
        5.,  // First coordinate maximum value
        -5., // Second coordinate minimum value
        5.   // Second coordinate maximum value
    };

    std::vector<double> m_spacing = {0.5, 0.5};

    std::vector<BoundaryCondition> m_boundaryConditions = {
        BoundaryCondition::DIRICHLET, // First coordinate lower boundary condition
        BoundaryCondition::DIRICHLET, // First coordinate upper boundary condition
        BoundaryCondition::DIRICHLET, // Second coordinate lower boundary condition
        BoundaryCondition::DIRICHLET  // Second coordinate upper boundary condition
    };
};

class Profile {
public:
    Profile();
    ~Profile();

    unsigned short const GetGridDimension() const;
    Geometry const GetGridGeometry() const;
    std::vector<BoundaryCondition> const& GetGridBoundaryConditions() const;
    std::vector<double> const& GetGridLimits() const;
    std::vector<double> const& GetSpacing() const;

    void SetGridDimension(unsigned short const dimension);
    void SetGridGeometry(Geometry const geometry);
    void SetGridBoundaryConditions(std::vector<BoundaryCondition> const& boundaryConditions);
    void SetGridLimits(std::vector<double> const& limits);
    void SetSpacing(std::vector<double> const& spacing);

private:
    std::unique_ptr<GridProfile> m_gridProfile;
};

}

#include "profile.hpp"
#include "profile_options.hpp"

#include <vector>

namespace MHD {

Profile::Profile() {
    m_gridDimension = Dimension::TWO;
    m_gridGeometry = Geometry::CARTESIAN;
    m_gridBounds = {-5.0, 5.0, -5.0, 5.0};
    m_gridSpacing = {0.5, 0.5};
    m_gridBoundaryConditions = {BoundaryCondition::DIRICHLET, BoundaryCondition::DIRICHLET, BoundaryCondition::DIRICHLET, BoundaryCondition::DIRICHLET};
}

Dimension const Profile::GetGridDimension() const {return m_gridDimension;}
Geometry const Profile::GetGridGeometry() const {return m_gridGeometry;}
std::vector<double> const& Profile::GetGridBounds() const {return m_gridBounds;}
std::vector<double> const& Profile::GetGridSpacing() const {return m_gridSpacing;}
std::vector<BoundaryCondition> const& Profile::GetGridBoundaryConditions() const {return m_gridBoundaryConditions;}

void Profile::SetGridDimension(Dimension const dimension) {m_gridDimension = dimension;}
void Profile::SetGridGeometry(Geometry const geometry) {m_gridGeometry = geometry;}
void Profile::SetGridBounds(std::vector<double> const& bounds) {m_gridBounds = bounds;}
void Profile::SetGridSpacing(std::vector<double> const& spacing) {m_gridSpacing = spacing;}
void Profile::SetGridBoundaryConditions(std::vector<BoundaryCondition> const& boundaryConditions) {m_gridBoundaryConditions = boundaryConditions;}

} // namespace MHD

#include "profile.hpp"

#include <memory>
#include <vector>

namespace MHD {

Profile::Profile() {
    m_gridProfile = std::make_unique<GridProfile>();
}

Profile::~Profile() {}

unsigned short const Profile::GetGridDimension() const {
    return m_gridProfile->m_dimension;
}

Geometry const Profile::GetGridGeometry() const {
    return m_gridProfile->m_geometry;
}

std::vector<BoundaryCondition> const& Profile::GetGridBoundaryConditions() const {
    return m_gridProfile->m_boundaryConditions;
}

std::vector<double> const& Profile::GetGridLimits() const {
    return m_gridProfile->m_limits;
}

std::vector<double> const& Profile::GetSpacing() const {
    return m_gridProfile->m_spacing;
}

void Profile::SetGridDimension(unsigned short const dimension) {
    m_gridProfile->m_dimension = dimension;
}

void Profile::SetGridGeometry(Geometry const geometry) {
    m_gridProfile->m_geometry = geometry;
}

void Profile::SetGridBoundaryConditions(std::vector<BoundaryCondition> const& boundaryConditions) {
    m_gridProfile->m_boundaryConditions = boundaryConditions;
}

void Profile::SetGridLimits(std::vector<double> const& limits) {
    m_gridProfile->m_limits = limits;
}

void Profile::SetSpacing(std::vector<double> const& spacing) {
    m_gridProfile->m_spacing = spacing;
}

}

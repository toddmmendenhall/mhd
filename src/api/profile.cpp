#include "profile.hpp"

#include <array>
#include <memory>

namespace MHD {

Profile::Profile() {
    m_gridProfile = std::make_unique<GridProfile>();
}

Profile::~Profile() {}

short const Profile::GetGridDimension() const {
    return m_gridProfile->m_dimension;
}

Geometry const Profile::GetGridGeometry() const {
    return m_gridProfile->m_geometry;
}

std::array<BoundaryCondition, 6> const& Profile::GetGridBoundaryConditions() const {
    return m_gridProfile->m_boundaryConditions;
}

std::array<double, 6> const& Profile::GetGridLimits() const {
    return m_gridProfile->m_limits;
}

void Profile::SetGridDimension(short const dimension) {
    m_gridProfile->m_dimension = dimension;
}

void Profile::SetGridGeometry(Geometry const geometry) {
    m_gridProfile->m_geometry = geometry;
}

void Profile::SetGridBoundaryConditions(std::array<BoundaryCondition, 6> const& boundaryConditions) {
    m_gridProfile->m_boundaryConditions = boundaryConditions;
}

void Profile::SetGridLimits(std::array<double, 6> const& limits) {
    m_gridProfile->m_limits = limits;
}

}

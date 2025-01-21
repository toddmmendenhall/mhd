#include "profile.hpp"
#include "profile_options.hpp"

#include <vector>

namespace MHD {

Profile::Profile() {
    m_gridDimension = Dimension::ONE;
    m_gridGeometry = Geometry::CARTESIAN;
    m_gridBounds = {-1.0, 1.0};
    m_gridSpacing = {0.2};
    m_gridBoundaryConditions = BoundaryCondition::DIRICHLET;
    m_spatialDerivativeMethod = SpatialDerivativeMethod::FINITE_DIFFERENCE;
    m_temporalIntegrationMethod = TemporalIntegrationMethod::FORWARD_EULER;
}

Dimension const Profile::GetGridDimension() const {
    return m_gridDimension;
}

Geometry const Profile::GetGridGeometry() const {
    return m_gridGeometry;
}

std::vector<double> const& Profile::GetGridBounds() const {
    return m_gridBounds;
}

std::vector<double> const& Profile::GetGridSpacing() const {
    return m_gridSpacing;
}

BoundaryCondition const& Profile::GetGridBoundaryConditions() const {
    return m_gridBoundaryConditions;
}

SpatialDerivativeMethod const Profile::GetSpatialDerivativeMethod() const {
    return m_spatialDerivativeMethod;
}

TemporalIntegrationMethod const Profile::GetTemporalIntegrationMethod() const {
    return m_temporalIntegrationMethod;
}

void Profile::SetGridDimension(Dimension const dimension) {
    m_gridDimension = dimension;
}

void Profile::SetGridGeometry(Geometry const geometry) {
    m_gridGeometry = geometry;
}

void Profile::SetGridBounds(std::vector<double> const& bounds) {
    m_gridBounds = bounds;
}

void Profile::SetGridSpacing(std::vector<double> const& spacing) {
    m_gridSpacing = spacing;
}

void Profile::SetGridBoundaryConditions(BoundaryCondition const& boundaryConditions) {
    m_gridBoundaryConditions = boundaryConditions;
}

void Profile::SetSpatialDerivativeMethod(SpatialDerivativeMethod const method) {
    m_spatialDerivativeMethod = method;
}

void Profile::SetTemporalIntegrationMethod(TemporalIntegrationMethod const method) {
    m_temporalIntegrationMethod = method;
}

} // namespace MHD

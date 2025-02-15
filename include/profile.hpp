#pragma once

#include <profile_options.hpp>

#include <vector>

namespace MHD {

struct Profile {
    // Grid options
    Dimension m_gridDimensionOption = Dimension::ONE;
    Geometry m_gridGeometryOption = Geometry::CARTESIAN;
    std::vector<double> m_gridBoundsOption = {0.0, 1.0, 0.0, 1.0, 0.0, 1.0};
    std::vector<double> m_gridSpacingsOption = {0.1, 0.1, 0.1};

    // Solver options
    BoundaryCondition m_gridBoundaryConditionOption = BoundaryCondition::DIRICHLET;
    SpatialDerivativeMethod m_spatialDerivativeOption = SpatialDerivativeMethod::FINITE_DIFFERENCE;
    TemporalIntegrationMethod m_temporalIntegrationOption = TemporalIntegrationMethod::FORWARD_EULER;
    FluxScheme m_fluxOption = FluxScheme::HLLC;
    ReconstructionOption m_reconstructionOption = ReconstructionOption::LINEAR;

    // Phenomenon options
    CompressibleOption m_compressibleOption = CompressibleOption::COMPRESSIBLE;
    ViscousOption m_viscousOption = ViscousOption::INVISCID;
};

} // namespace MHD

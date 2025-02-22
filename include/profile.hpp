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
    BoundaryConditionOption m_boundaryConditionOption = BoundaryConditionOption::OUTFLOW;
    SpatialDerivativeMethod m_spatialDerivativeOption = SpatialDerivativeMethod::FINITE_DIFFERENCE;
    TemporalIntegrationMethod m_temporalIntegrationOption = TemporalIntegrationMethod::FORWARD_EULER;
    FluxScheme m_fluxOption = FluxScheme::LOG;
    ReconstructionOption m_reconstructionOption = ReconstructionOption::CONSTANT;

    // Phenomenon options
    CompressibleOption m_compressibleOption = CompressibleOption::COMPRESSIBLE;
    ViscousOption m_viscousOption = ViscousOption::INVISCID;
    HydroOption m_hydroOption = HydroOption::PURE_HYDRO;

    // Generic options
    OutputDataOption m_outputDataOption = OutputDataOption::NO;
};

} // namespace MHD

#pragma once

#include <profile_options.hpp>

#include <vector>

namespace MHD {

struct Profile {
    // Grid options
    Dimension m_gridDimensionOption = Dimension::ONE;
    std::vector<double> m_gridBoundsOption = {0.0, 20.0, 0.0, 1.0, 0.0, 1.0};
    std::vector<double> m_gridSpacingsOption = {0.1, 0.1, 0.1};

    // Solver options
    BoundaryConditionOption m_boundaryConditionOption = BoundaryConditionOption::REFLECTIVE;
    ReconstructionOption m_reconstructionOption = ReconstructionOption::MUSCL;
    FluxScheme m_fluxOption = FluxScheme::KT;
    TemporalIntegrationMethod m_temporalIntegrationOption = TemporalIntegrationMethod::FORWARD_EULER;

    // Phenomenon options
    CompressibleOption m_compressibleOption = CompressibleOption::COMPRESSIBLE;

    // Generic options
    OutputDataOption m_outputDataOption = OutputDataOption::NO;
};

} // namespace MHD

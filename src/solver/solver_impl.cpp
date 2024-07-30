#include "grid.hpp"
#include "profile.hpp"
#include "solver_impl.hpp"

#include <memory>

namespace MHD {

SolverImpl::SolverImpl(std::unique_ptr<Profile> const& profile) {
    profile->GetSpatialDerivativeMethod();
    profile->GetTemporalIntegrationMethod();
}

SolverImpl::~SolverImpl() {}

void SolverImpl::SolveTimeStep(std::unique_ptr<Grid> const& grid) {
    grid->GetGridImpl()->;
}
    
} // namespace MHD

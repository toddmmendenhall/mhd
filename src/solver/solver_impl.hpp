#pragma once

#include "grid.hpp"
#include "profile.hpp"
#include "solver_impl.hpp"

#include <memory>

namespace MHD {

class SolverImpl {
public:
    SolverImpl(std::unique_ptr<Profile> const& profile);
    ~SolverImpl();

    void SolveTimeStep(std::unique_ptr<Grid> const& grid);
};

} // namespace MHD

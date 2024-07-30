#pragma once

#include "grid.hpp"
#include "profile.hpp"

#include <memory>

namespace MHD {

class SolverImpl;

class Solver {
public:
    Solver(std::unique_ptr<Profile> const& profile);
    ~Solver();

    void SolveTimeStep(std::unique_ptr<Grid> const& grid);

private:
    std::unique_ptr<SolverImpl> m_solverImpl;
};
    
} // namespace MHD

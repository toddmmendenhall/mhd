#pragma once

#include "grid.hpp"
#include "profile.hpp"
#include "solver.hpp"

#include <memory>

namespace MHD {

class Calc {
public:
    Calc();

    void Run();

    std::unique_ptr<Profile> m_profile;
    std::unique_ptr<Grid> m_grid;
    std::unique_ptr<Solver> m_solver;
    
private:
    void SetupCalc();

};

} // namespace MHD

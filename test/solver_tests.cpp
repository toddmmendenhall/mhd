#include <profile.hpp>
#include <profile_options.hpp>
#include <solver.hpp>
#include <cartesian/1d.hpp>

#include "gtest/gtest.h"

#include <iostream>

TEST(SolverTests, SpatialDerivative) {
    MHD::Profile profile;
    MHD::Cartesian1DGrid grid(profile);
    MHD::Solver solver(profile);
    solver.SolveTimeStep(grid);
}

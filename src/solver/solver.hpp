#pragma once

#include <memory>

namespace MHD {

class Profile;
class Grid;
class State;

class Solver {
public:
    Solver(Profile const& profile);
    ~Solver();

    void SolveTimeStep(Grid const& grid);

private:
    std::unique_ptr<State> m_state;
};

} // namespace MHD

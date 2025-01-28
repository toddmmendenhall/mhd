#pragma once

#include <memory>

namespace MHD {

class IGrid;
class Profile;
class Solver;

class Calc {
public:
    Calc(Profile const& profile);
    ~Calc() = default;

    void Run();
    
private:
    Profile const& m_profile;
    std::unique_ptr<IGrid> m_grid;
    std::unique_ptr<Solver> m_solver;
};

} // namespace MHD

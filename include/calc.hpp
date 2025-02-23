#pragma once

#include <memory>

namespace MHD {

class ExecutionController;
class IGrid;
class ISolver;
class Profile;
class VariableStore;

enum class InitialCondition {
    SOD_SHOCK_TUBE = 0,
};

class Calc {
public:
    Calc(Profile const& profile);
    ~Calc();

    void SetInitialCondition(InitialCondition ic);
    void Run();
    
private:
    void SetSodShockTube();

    void WriteData(VariableStore const& varStore);

    std::unique_ptr<ExecutionController> m_executionController;
    std::unique_ptr<IGrid> m_grid;
    Profile const& m_profile;
    std::unique_ptr<ISolver> m_solver;
    std::unique_ptr<VariableStore> m_variableStore;
    double const m_duration = 1e-1;
    double const tStep = 1e-6;
    double m_currentTime = 0.0;
    std::size_t m_currentStep = 0;
    std::size_t m_currentOutput = 0;
};

} // namespace MHD

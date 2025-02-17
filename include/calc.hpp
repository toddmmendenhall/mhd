#pragma once

#include <memory>

namespace MHD {

class ExecutionController;
class IGrid;
class IIntegrator;
class ISolver;
class Profile;
class VariableStore;

class Calc {
public:
    Calc(Profile const& profile);
    ~Calc();

    void SetInitialConditions();
    void SetSodShockTube();
    void Run();
    
private:
    void WriteData(VariableStore const& varStore);

    std::unique_ptr<ExecutionController> m_executionController;
    std::unique_ptr<IGrid> m_grid;
    Profile const& m_profile;
    std::unique_ptr<IIntegrator> m_integrator;
    std::unique_ptr<ISolver> m_solver;
    std::unique_ptr<VariableStore> m_variableStore;
    double const m_duration = 0.001;
    double m_currentTime = 0.0;
    std::size_t m_currentStep = 0;
};

} // namespace MHD

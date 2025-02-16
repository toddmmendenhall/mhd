#pragma once

#include <memory>

namespace MHD {

class ExecutionController;
class IGrid;
class ISolver;
class Profile;
class VariableStore;

class Calc {
public:
    Calc(Profile const& profile);
    ~Calc();

    void Run();
    
private:
    void WriteData(VariableStore const& varStore);

    std::unique_ptr<ExecutionController> m_executionController;
    std::unique_ptr<IGrid> m_grid;
    Profile const& m_profile;
    std::unique_ptr<ISolver> m_solver;
    std::unique_ptr<VariableStore> m_variableStore;
};

} // namespace MHD

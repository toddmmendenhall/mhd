#pragma once

#include <error.hpp>

#include <memory>

namespace MHD {

class ExecutionController;
class IFluxScheme;
class IReconstruction;
class Profile;
class VariableStore;

class ISolver {
public:
    virtual ~ISolver() = default;
    virtual Error PerformTimeStep() = 0;
};

class Solver : public ISolver {
public:
    Solver(Profile const& profile, ExecutionController const& execCtrl,
           VariableStore& varStore);
    
    Error PerformTimeStep();

    void ComputePrimitivesFromConserved();

    void UpdateConservedFromPrimitives();

    void ComputeFluxes();

    void ReconstructVariables();

private:
    std::unique_ptr<IFluxScheme> m_fluxScheme;
    std::unique_ptr<IReconstruction> m_reconstruction;
    ExecutionController const& m_execCtrl;
    VariableStore& m_varStore;
};

std::unique_ptr<ISolver> solverFactory(Profile const& profile, ExecutionController const& execCtrl,
                                       VariableStore& varStore);

} // namespace MHD

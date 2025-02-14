#pragma once

#include <error.hpp>

#include <memory>

namespace MHD {

class ExecutionController;
class Profile;
class VariableStore;
class IFluxScheme;
class IReconstruction;

class ISolver {
public:
    virtual ~ISolver() = default;
    virtual Error PerformTimeStep() = 0;
};

class Solver : public ISolver {
public:
    Solver(Profile const& profile);
    
    Error PerformTimeStep();

    void ComputePrimitivesFromConserved();

    void UpdateConservedFromPrimitives();

    void ComputeFluxes();

    void ReconstructVariables();

private:
    std::unique_ptr<VariableStore> varStore;
    std::unique_ptr<ExecutionController> execCtrl;
    std::unique_ptr<IFluxScheme> m_fluxScheme;
    std::unique_ptr<IReconstruction> m_reconstruction;
};

std::unique_ptr<ISolver> solverFactory(Profile const& profile);

} // namespace MHD

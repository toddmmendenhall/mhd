#pragma once

#include <memory>

namespace MHD {

class ExecutionController;
class Profile;
class VariableStore;
class IFluxScheme;
class IReconstruction;

class Solver {
public:
    Solver(Profile const& profile);
    ~Solver();

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

} // namespace MHD

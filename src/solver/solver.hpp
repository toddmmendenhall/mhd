#pragma once

#include <memory>

namespace MHD {

class ExecutionController;
class Profile;
class VariableStore;
class IFluxScheme;
class IGrid;

class Solver {
public:
    Solver(Profile const& profile, IGrid const& grid);
    ~Solver();

    void ComputePrimitivesFromConserved();

    void UpdateConservedFromPrimitives();

    void ComputeFluxes();

private:
    std::unique_ptr<VariableStore> varStore;
    std::unique_ptr<ExecutionController> execCtrl;
    IGrid const& m_grid;
    std::unique_ptr<IFluxScheme> m_fluxScheme;
};

} // namespace MHD

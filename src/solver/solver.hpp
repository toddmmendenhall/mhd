#pragma once

#include <error.hpp>

#include <memory>

namespace MHD {

class ElectricFieldCalculator;
class ExecutionController;
class IBoundaryCondition;
class IFluxScheme;
class IGrid;
class IReconstruction;
class MagneticFieldCalculator;
class Profile;
class VariableStore;

class ISolver {
public:
    virtual ~ISolver() = default;
    virtual Error PerformTimeStep(IGrid const& grid) = 0;
};

class Solver : public ISolver {
public:
    Solver(Profile const& profile, ExecutionController const& execCtrl,
           VariableStore& varStore);
    
    Error PerformTimeStep(IGrid const& grid);

    void ComputePrimitivesFromConserved();

    void UpdateConservedFromPrimitives();

    void ComputeFluxes(IGrid const& grid);

    void ReconstructVariables();

    void ComputeElectricFields();

    void ComputeMagneticFields();

    void ApplyBoundaryConditions();

private:
    std::unique_ptr<MagneticFieldCalculator> m_bFieldCalc;
    std::unique_ptr<ElectricFieldCalculator> m_eFieldCalc;
    std::unique_ptr<IBoundaryCondition> m_boundCon;
    std::unique_ptr<IFluxScheme> m_fluxScheme;
    std::unique_ptr<IReconstruction> m_reconstruction;
    ExecutionController const& m_execCtrl;
    VariableStore& m_varStore;
};

std::unique_ptr<ISolver> solverFactory(Profile const& profile, ExecutionController const& execCtrl,
                                       VariableStore& varStore);

} // namespace MHD

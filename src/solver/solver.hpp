#pragma once

#include <memory>

namespace MHD {

class BoundaryConditionContext;
class ElectricFieldCalculator;
class ExecutionController;
struct FluxContext;
class IBoundaryCondition;
class IFluxScheme;
class IGrid;
class IIntegrator;
struct IntegrationContext;
class IReconstruction;
class MagneticFieldCalculator;
class Profile;
struct ReconstructionContext;
struct ResidualContext;
class Residual;
struct TransportContext;
class VariableStore;

class ISolver {
public:
    virtual ~ISolver() = default;
    virtual void SetupConservedState() = 0;
    virtual void PerformTimeStep() = 0;
    virtual double const TimeStep() const = 0;
};

class Solver : public ISolver {
public:
    Solver(Profile const& profile, ExecutionController const& execCtrl,
           VariableStore& varStore, IGrid const& grid);
    
    void SetupConservedState();
    
    void PerformTimeStep();

    void ComputePrimitivesFromConserved();

    void ComputeFluxes();

    void ReconstructVariables();

    void ApplyBoundaryConditions();

    void ComputeResiduals();

    void CalculateTimeStep();

    double const TimeStep() const { return timeStep; }

private:
    double cfl = 0.5;
    double timeStep = 1e-5;
    std::unique_ptr<FluxContext> m_fluxContext;
    std::unique_ptr<IBoundaryCondition> m_boundCon;
    std::unique_ptr<IIntegrator> m_integrator;
    std::unique_ptr<IFluxScheme> m_fluxScheme;
    std::unique_ptr<ReconstructionContext> m_reconstructionContext;
    std::unique_ptr<IReconstruction> m_reconstruction;
    std::unique_ptr<ResidualContext> m_residualContext;
    std::unique_ptr<Residual> m_residual;
    std::unique_ptr<IntegrationContext> m_integrationContext;
    std::unique_ptr<BoundaryConditionContext> m_boundaryConditionContext;
    ExecutionController const& m_execCtrl;
    IGrid const& m_grid;
    VariableStore& m_varStore;
};

std::unique_ptr<ISolver> solverFactory(Profile const& profile, ExecutionController const& execCtrl,
                                       VariableStore& varStore, IGrid const& grid);

} // namespace MHD

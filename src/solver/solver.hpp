#pragma once

#include <memory>

namespace MHD {

class ExecutionController;
class IBoundaryCondition;
class IFlux;
class IGrid;
class IIntegrator;
class IReconstruction;
class Profile;
class Residual;
class VariableStore;

class ISolver {
public:
    virtual ~ISolver() = default;
    virtual void PerformTimeStep() = 0;
    virtual double const TimeStep() const = 0;
    virtual void PrimFromCons() = 0;
};

class Solver : public ISolver {
public:
    Solver(Profile const& profile, ExecutionController const& execCtrl,
           VariableStore& varStore, IGrid const& grid);
    
    void ConsFromPrim();
    
    void PerformTimeStep();

    void PrimFromCons();

    void CalculateTimeStep();

    double const TimeStep() const { return timeStep; }

private:
    double cfl = 0.4;
    double timeStep = 1e-5;
    std::unique_ptr<IBoundaryCondition> m_boundCon;
    std::unique_ptr<IIntegrator> m_integrator;
    std::unique_ptr<IFlux> m_flux;
    std::unique_ptr<IReconstruction> m_reconstruction;
    std::unique_ptr<Residual> m_residual;
    ExecutionController const& m_execCtrl;
    IGrid const& m_grid;
    VariableStore& m_varStore;
};

std::unique_ptr<ISolver> solverFactory(Profile const& profile, ExecutionController const& execCtrl,
                                       VariableStore& varStore, IGrid const& grid);

} // namespace MHD

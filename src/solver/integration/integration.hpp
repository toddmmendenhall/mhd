#pragma once

#include <execution_controller.hpp>
#include <variable_store.hpp>
#include <residual.hpp>

#include <memory>

namespace MHD {

struct IntegrationContext {
    IntegrationContext(ResidualContext const& rc, VariableStore& vs, double const& tStep) : 
        tStep(tStep), cellIdxs(rc.cellIdxs), rhoRes(rc.rhoRes), rhoURes(rc.rhoURes), rhoVRes(rc.rhoVRes),
        rhoWRes(rc.rhoWRes), rhoERes(rc.rhoERes), rho(vs.rho), rhoU(vs.rhoU), rhoV(vs.rhoV), rhoW(vs.rhoW),
        rhoE(vs.rhoE) {}

    double const& tStep;
    std::vector<std::size_t> const& cellIdxs;

    // Cell-centered residuals
    std::vector<double> const& rhoRes;     // mass density residual
    std::vector<double> const& rhoURes;    // x momentum density residual
    std::vector<double> const& rhoVRes;    // y momentum density residual
    std::vector<double> const& rhoWRes;    // z momentum density residual
    std::vector<double> const& rhoERes;    // total energy density residual

    // Cell-centered states
    std::vector<double>& rho;
    std::vector<double>& rhoU;
    std::vector<double>& rhoV;
    std::vector<double>& rhoW;
    std::vector<double>& rhoE;
};

class IIntegrator {
public:
    virtual ~IIntegrator() = default;
    virtual void Compute(ExecutionController const& execCtrl) = 0;
    IntegrationContext const& GetContext() const { return *m_context; }

protected:
    std::unique_ptr<IntegrationContext> m_context;
};

struct ForwardEulerKernel {
    ForwardEulerKernel(IntegrationContext const& context) : m_context(context) {}

    void operator()(std::size_t const i) {
        m_context.rho[i] += m_context.tStep * m_context.rhoRes[i];
        m_context.rhoU[i] += m_context.tStep * m_context.rhoURes[i];
        m_context.rhoV[i] += m_context.tStep * m_context.rhoVRes[i];
        m_context.rhoW[i] += m_context.tStep * m_context.rhoWRes[i];
        m_context.rhoE[i] += m_context.tStep * m_context.rhoERes[i];
    }

    IntegrationContext const& m_context;
};

class ForwardEuler : public IIntegrator {
public:
    ForwardEuler(ResidualContext const& rc, VariableStore& vs, double const& tStep) {
        m_context = std::make_unique<IntegrationContext>(rc, vs, tStep);
    }

    void Compute(ExecutionController const& execCtrl) {
        ForwardEulerKernel kern(*m_context);
        execCtrl.LaunchKernel(kern, m_context->cellIdxs);
    }
};

std::unique_ptr<IIntegrator> integratorFactory(ResidualContext const& rc, VariableStore& vs, double const& tStep) {
    return std::make_unique<ForwardEuler>(rc, vs, tStep);
}

} // namespace MHD
#pragma once

#include <execution_controller.hpp>
#include <variable_store.hpp>
#include <residual.hpp>

#include <memory>

namespace MHD {

struct IntegrationContext {
    IntegrationContext(ResidualContext const& rc, VariableStore& vs, double const& tStep) : 
        tStep(tStep), numCells(rc.numCells), rhoRes(rc.rhoRes), rhoURes(rc.rhoURes), rhoVRes(rc.rhoVRes),
        rhoWRes(rc.rhoWRes), rhoERes(rc.rhoERes), bxRes(rc.bxRes), byRes(rc.byRes), bzRes(rc.bzRes),
        rho(vs.rho), rhoU(vs.rhoU), rhoV(vs.rhoV), rhoW(vs.rhoW),
        rhoE(vs.rhoE), bx(vs.bx), by(vs.by), bz(vs.bz) {}

    double const& tStep;
    std::size_t const numCells;

    // Cell-centered residuals
    std::vector<double> const& rhoRes;     // mass density residual
    std::vector<double> const& rhoURes;    // x momentum density residual
    std::vector<double> const& rhoVRes;    // y momentum density residual
    std::vector<double> const& rhoWRes;    // z momentum density residual
    std::vector<double> const& rhoERes;    // total energy density residual
    std::vector<double> const& bxRes;       // x magnetic field residual
    std::vector<double> const& byRes;       // y magnetic field residual
    std::vector<double> const& bzRes;       // z magnetic field residual

    // Cell-centered states
    std::vector<double>& rho;
    std::vector<double>& rhoU;
    std::vector<double>& rhoV;
    std::vector<double>& rhoW;
    std::vector<double>& rhoE;
    std::vector<double>& bx;
    std::vector<double>& by;
    std::vector<double>& bz;
};

class IIntegrator {
public:
    virtual ~IIntegrator() = default;
    virtual void Integrate(ExecutionController const& execCtrl) = 0;
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
        m_context.bx[i] += m_context.tStep * m_context.bxRes[i];
        m_context.by[i] += m_context.tStep * m_context.byRes[i];
        m_context.bz[i] += m_context.tStep * m_context.bzRes[i];
    }

    IntegrationContext const& m_context;
};

class ForwardEuler : public IIntegrator {
public:
    ForwardEuler(ResidualContext const& rc, VariableStore& vs, double const& tStep) {
        m_context = std::make_unique<IntegrationContext>(rc, vs, tStep);
    }

    void Integrate(ExecutionController const& execCtrl) {
        ForwardEulerKernel kern(*m_context);
        execCtrl.LaunchKernel(kern, m_context->numCells);
    }
};

std::unique_ptr<IIntegrator> integratorFactory(ResidualContext const& rc, VariableStore& vs, double const& tStep) {
    return std::make_unique<ForwardEuler>(rc, vs, tStep);
}

} // namespace MHD
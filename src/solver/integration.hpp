#pragma once

#include <context.hpp>
#include <execution_controller.hpp>
#include <variable_store.hpp>

#include <memory>

namespace MHD {

class IIntegrator {
public:
    virtual ~IIntegrator() = default;
    virtual void Compute(ExecutionController const& execCtrl, IntegrationContext const& context) = 0;
};

struct ForwardEulerKernel {
    ForwardEulerKernel(IntegrationContext const& context) : m_context(context) {}

    void operator()(std::size_t i) {
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
    void Compute(ExecutionController const& execCtrl, IntegrationContext const& context) {
        ForwardEulerKernel kern(context);
        execCtrl.LaunchKernel(kern, context.numCells);
    }
};

std::unique_ptr<IIntegrator> integratorFactory() {
    return std::make_unique<ForwardEuler>();
}

} // namespace MHD
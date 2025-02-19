#pragma once

#include <context.hpp>
#include <execution_controller.hpp>
#include <variable_store.hpp>

#include <memory>

namespace MHD {

class IIntegrator {
public:
    virtual ~IIntegrator() = default;
    virtual void Solve(ExecutionController const& execCtrl, IntegrationContext const& context, VariableStore& varStore) = 0;
};

struct ForwardEulerKernel {
    ForwardEulerKernel(IntegrationContext const& context, VariableStore& varStore)
        : m_context(context), m_varStore(varStore) {}

    void operator()(std::size_t i) {
        m_varStore.rho[i] += m_context.timeStep * m_context.rhoRes[i];
        m_varStore.rhoU[i] += m_context.timeStep * m_context.rhoURes[i];
        m_varStore.rhoV[i] += m_context.timeStep * m_context.rhoVRes[i];
        m_varStore.rhoW[i] += m_context.timeStep * m_context.rhoWRes[i];
        m_varStore.rhoE[i] += m_context.timeStep * m_context.rhoERes[i];
    }

    IntegrationContext const& m_context;
    VariableStore& m_varStore;
};

class ForwardEuler : public IIntegrator {
public:
    void Solve(ExecutionController const& execCtrl, IntegrationContext const& context, VariableStore& varStore) override {
        ForwardEulerKernel kern(context, varStore);
        execCtrl.LaunchKernel(kern, varStore.numCells);
    }
};

std::unique_ptr<IIntegrator> integratorFactory() {
    return std::make_unique<ForwardEuler>();
}

} // namespace MHD
#pragma once

#include <execution_controller.hpp>
#include <variable_store.hpp>

#include <memory>

namespace MHD {

class IIntegrator {
public:
    virtual ~IIntegrator() = default;
    virtual void Solve(ExecutionController const& execCtrl, VariableStore& varStore) = 0;
};

struct ForwardEulerKernel {
    ForwardEulerKernel(VariableStore& varStore) : m_varStore(varStore) {}

    void operator()(std::size_t idx) {
        double const cfl = 0.5;
        double dRho = 1.0;
        if (idx == 0) {
            dRho = m_varStore.m_rho[idx] - m_varStore.m_rho[10];
        } else if (idx == 9) {
            dRho = m_varStore.m_rho[11] - m_varStore.m_rho[idx];
        } else {
            dRho = m_varStore.m_rho[idx] - m_varStore.m_rho[idx - 1];
        }
        m_varStore.m_rho[idx] -= cfl * dRho;
    }
    VariableStore& m_varStore;
};

class ForwardEuler : public IIntegrator {
public:
    void Solve(ExecutionController const& execCtrl, VariableStore& varStore) {
        ForwardEulerKernel kern(varStore);
        execCtrl.LaunchKernel(kern, 10);
    }

};

std::unique_ptr<IIntegrator> integratorFactory() {
    return std::make_unique<ForwardEuler>();
}

} // namespace MHD
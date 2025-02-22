#pragma once

#include <context.hpp>
#include <execution_controller.hpp>

namespace MHD {

struct TransportKernel {
public:
    TransportKernel(ResidualContext& context) : m_context(context) {}

    void operator()(std::size_t i) {
        // Get the left and right face indices for this cell
        std::size_t iLeft = m_context.cellToFaceIndices[i][0];
        std::size_t iRight = m_context.cellToFaceIndices[i][1];
        double c = -1.0 / m_context.cellSize[0];

        m_context.rhoRes[i] = c * (m_context.rhoFlux[iRight] - m_context.rhoFlux[iLeft]);
        m_context.rhoURes[i] = c * (m_context.rhoUFlux[iRight] - m_context.rhoUFlux[iLeft]);
        m_context.rhoVRes[i] = c * (m_context.rhoVFlux[iRight] - m_context.rhoVFlux[iLeft]);
        m_context.rhoWRes[i] = c * (m_context.rhoWFlux[iRight] - m_context.rhoWFlux[iLeft]);
        m_context.rhoERes[i] = c * (m_context.rhoEFlux[iRight] - m_context.rhoEFlux[iLeft]);
    }

    ResidualContext& m_context;
};

class Residual {
public:
    ~Residual() = default;
    void Compute(ExecutionController const& execCtrl, ResidualContext& context) {
        TransportKernel kernel(context);
        execCtrl.LaunchKernel(kernel, context.numCells);
    }
};

} // namespace MHD
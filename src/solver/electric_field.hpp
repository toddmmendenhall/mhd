#pragma once

#include <context.hpp>
#include <execution_controller.hpp>

namespace MHD {

struct ElectricFieldKernel {
    ElectricFieldKernel(ElectricFieldContext& context) : m_context(context) {}
    inline void operator()(std::size_t const edgeIdx) {
        m_context.eZ[edgeIdx] = 0.25 * ( -m_context.fluxBY[m_context.edgeToFaceIdx[edgeIdx][0]]
                                          - m_context.fluxBY[m_context.edgeToFaceIdx[edgeIdx][1]]
                                          + m_context.fluxBX[m_context.edgeToFaceIdx[edgeIdx][2]]
                                          + m_context.fluxBX[m_context.edgeToFaceIdx[edgeIdx][3]]);
    }
    ElectricFieldContext& m_context;
};

struct MagneticFieldKernel {
    MagneticFieldKernel(MagneticFieldContext& context) : m_context(context) {}
    inline void operator()(std::size_t const faceIdx) {
        m_context.bX[faceIdx] -= m_context.timeStep / m_context.cellSize[1] * (
                                   m_context.eZ[m_context.faceToEdgeIdx[faceIdx][0]] -
                                   m_context.eZ[m_context.faceToEdgeIdx[faceIdx][1]]);
        m_context.bY[faceIdx] -= m_context.timeStep / m_context.cellSize[0] * (
                                   -m_context.eZ[m_context.faceToEdgeIdx[faceIdx][0]] +
                                   m_context.eZ[m_context.faceToEdgeIdx[faceIdx][1]]);
    }
    MagneticFieldContext& m_context;
};

class ElectricFieldCalculator {
public:
    ElectricFieldCalculator() = default;
    ~ElectricFieldCalculator() = default;

    void Compute(ExecutionController const& execCtrl, ElectricFieldContext& context) {
        ElectricFieldKernel kern(context);
        execCtrl.LaunchKernel(kern, context.edgeToFaceIdx.size());
    }
};

class MagneticFieldCalculator {
    public:
        MagneticFieldCalculator() = default;
        ~MagneticFieldCalculator() = default;
    
        void Compute(ExecutionController const& execCtrl, MagneticFieldContext& context) {
            MagneticFieldKernel kern(context);
            execCtrl.LaunchKernel(kern, context.faceToEdgeIdx.size());
        }
    };

} // namespace MHD
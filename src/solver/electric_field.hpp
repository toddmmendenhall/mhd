#pragma once

#include <context.hpp>
#include <execution_controller.hpp>

namespace MHD {

struct ElectricFieldKernel {
    ElectricFieldKernel(ElectricFieldContext& context) : m_context(context) {}
    inline void operator()(std::size_t const edgeIdx) {
        m_context.m_eZ[edgeIdx] = 0.25 * ( -m_context.m_fluxBY[m_context.m_edgeToFaceIdx[edgeIdx][0]]
                                          - m_context.m_fluxBY[m_context.m_edgeToFaceIdx[edgeIdx][1]]
                                          + m_context.m_fluxBX[m_context.m_edgeToFaceIdx[edgeIdx][2]]
                                          + m_context.m_fluxBX[m_context.m_edgeToFaceIdx[edgeIdx][3]]);
    }
    ElectricFieldContext& m_context;
};

struct MagneticFieldKernel {
    MagneticFieldKernel(MagneticFieldContext& context) : m_context(context) {}
    inline void operator()(std::size_t const faceIdx) {
        m_context.m_bX[faceIdx] -= m_context.m_timeStep / m_context.m_cellSize[1] * (
                                   m_context.m_eZ[m_context.m_faceToEdgeIdx[faceIdx][0]] -
                                   m_context.m_eZ[m_context.m_faceToEdgeIdx[faceIdx][1]]);
        m_context.m_bY[faceIdx] -= m_context.m_timeStep / m_context.m_cellSize[0] * (
                                   -m_context.m_eZ[m_context.m_faceToEdgeIdx[faceIdx][0]] +
                                   m_context.m_eZ[m_context.m_faceToEdgeIdx[faceIdx][1]]);
    }
    MagneticFieldContext& m_context;
};

class ElectricFieldCalculator {
public:
    ElectricFieldCalculator() = default;
    ~ElectricFieldCalculator() = default;

    void Compute(ExecutionController const& execCtrl, ElectricFieldContext& context) {
        ElectricFieldKernel kern(context);
        execCtrl.LaunchKernel(kern, context.m_edgeToFaceIdx.size());
    }
};

class MagneticFieldCalculator {
    public:
        MagneticFieldCalculator() = default;
        ~MagneticFieldCalculator() = default;
    
        void Compute(ExecutionController const& execCtrl, MagneticFieldContext& context) {
            MagneticFieldKernel kern(context);
            execCtrl.LaunchKernel(kern, context.m_faceToEdgeIdx.size());
        }
    };

} // namespace MHD
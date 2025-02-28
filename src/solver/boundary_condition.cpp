#include <boundary_condition.hpp>
#include <context.hpp>
#include <execution_controller.hpp>
#include <profile.hpp>

#include <memory>

namespace {
 
struct OutflowBoundaryConditionKernel {
    OutflowBoundaryConditionKernel(MHD::BoundaryConditionContext& context) :
        m_context(context) {}

    void operator()(std::size_t i) {
        std::size_t const iBnd = m_context.boundaryFaceToBoundaryCellIndices[i];
        std::size_t const iInt = m_context.boundaryFaceToInteriorCellIndices[i];

        // Normal component of the velocity inside the boundary
        auto uBar = m_context.u[iInt] * m_context.faceNormalX[i] +
                    m_context.v[iInt] * m_context.faceNormalY[i] +
                    m_context.w[iInt] * m_context.faceNormalZ[i];

        // Neumann condition on the normal velocity
        m_context.u[iBnd] = uBar * m_context.faceNormalX[i];
        m_context.v[iBnd] = uBar * m_context.faceNormalY[i];
        m_context.w[iBnd] = uBar * m_context.faceNormalZ[i];

        // Dirichlet condition on the pressure
        m_context.rho[iBnd] = m_context.rho[iInt];
        m_context.p[iBnd] = m_context.p[iInt];
    }
    MHD::BoundaryConditionContext& m_context;
};

struct ReflectiveBoundaryConditionKernel {
    ReflectiveBoundaryConditionKernel(MHD::BoundaryConditionContext& context) :
        m_context(context) {}

    void operator()(std::size_t i) {
        std::size_t const iGhost = m_context.boundaryFaceToBoundaryCellIndices[i];
        std::size_t const iBoundary = m_context.boundaryFaceToInteriorCellIndices[i];

        // Copy scalars to ghost cells
        m_context.rho[iGhost] = m_context.rho[iBoundary];
        m_context.rho[iGhost+1] = m_context.rho[iBoundary];
        m_context.p[iGhost] = m_context.p[iBoundary];
        m_context.p[iGhost+1] = m_context.p[iBoundary];

        // Reflect vectors in ghost cells
        m_context.u[iGhost] = -m_context.u[iBoundary];
        m_context.v[iGhost] = -m_context.v[iBoundary];
        m_context.w[iGhost] = -m_context.w[iBoundary];
        m_context.u[iGhost+1] = -m_context.u[iBoundary];
        m_context.v[iGhost+1] = -m_context.v[iBoundary];
        m_context.w[iGhost+1] = -m_context.w[iBoundary];
    }

    MHD::BoundaryConditionContext& m_context;
};

} // namespace

namespace MHD {

class OutflowBoundaryCondition : public IBoundaryCondition {
public:
    OutflowBoundaryCondition() = default;

    void Compute(ExecutionController const& execCtrl, BoundaryConditionContext& context) {
        OutflowBoundaryConditionKernel kern(context);
        execCtrl.LaunchKernel(kern, context.boundaryFaceToBoundaryCellIndices.size());
    }
};

class ReflectiveBoundaryCondition : public IBoundaryCondition {
    public:
        ReflectiveBoundaryCondition() = default;
    
        void Compute(ExecutionController const& execCtrl, BoundaryConditionContext& context) {
            ReflectiveBoundaryConditionKernel kern(context);
            execCtrl.LaunchKernel(kern, context.boundaryFaceToBoundaryCellIndices.size());
        }
    };

std::unique_ptr<IBoundaryCondition> boundaryConditionFactory(Profile const& profile) {
    if (BoundaryConditionOption::REFLECTIVE == profile.m_boundaryConditionOption) {
        return std::make_unique<ReflectiveBoundaryCondition>();
    }
    if (BoundaryConditionOption::OUTFLOW == profile.m_boundaryConditionOption) {
        return std::make_unique<OutflowBoundaryCondition>();
    }
    return nullptr;
}

} // namespace MHD
#include <boundary_condition.hpp>
#include <context.hpp>
#include <execution_controller.hpp>
#include <profile.hpp>

#include <memory>

namespace {
 
struct OutflowBoundaryConditionKernel {
    OutflowBoundaryConditionKernel(MHD::BoundaryConditionContext& context) :
        m_context(context) {}

    void operator()(std::size_t idx) {
        // Normal component of the velocity inside the boundary
        auto uBar = m_context.uIn[idx] * m_context.faceNormalX[idx] +
                    m_context.vIn[idx] * m_context.faceNormalY[idx] +
                    m_context.wIn[idx] * m_context.faceNormalZ[idx];

        // Neumann condition on the normal velocity
        m_context.u[idx] = uBar * m_context.faceNormalX[idx];
        m_context.v[idx] = uBar * m_context.faceNormalY[idx];
        m_context.w[idx] = uBar * m_context.faceNormalZ[idx];

        // Dirichlet condition on the pressure
        m_context.p[idx] = m_context.pOut;
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
        execCtrl.LaunchKernel(kern, context.m_boundaryNodeIndices.size());
    }
};

std::unique_ptr<IBoundaryCondition> boundaryConditionFactory(Profile const& profile) {
    if (BoundaryConditionOption::OUTFLOW == profile.m_boundaryConditionOption) {
        return std::make_unique<OutflowBoundaryCondition>();
    }
    return nullptr;
}

} // namespace MHD
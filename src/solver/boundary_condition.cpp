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
        std::size_t const iBnd = m_context.boundaryToInteriorNodeIndices[i].first;
        std::size_t const iInt = m_context.boundaryToInteriorNodeIndices[i].second;
        std::size_t const iFace = m_context.boundaryToFaceIndices[i].second;

        // Normal component of the velocity inside the boundary
        auto uBar = m_context.u[iInt] * m_context.faceNormalX[iFace] +
                    m_context.v[iInt] * m_context.faceNormalY[iFace] +
                    m_context.w[iInt] * m_context.faceNormalZ[iFace];

        // Neumann condition on the normal velocity
        m_context.u[iBnd] = uBar * m_context.faceNormalX[iFace];
        m_context.v[iBnd] = uBar * m_context.faceNormalY[iFace];
        m_context.w[iBnd] = uBar * m_context.faceNormalZ[iFace];

        // Dirichlet condition on the pressure
        m_context.p[iBnd] = m_context.p[iInt];
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
        execCtrl.LaunchKernel(kern, context.boundaryToInteriorNodeIndices.size());
    }
};

std::unique_ptr<IBoundaryCondition> boundaryConditionFactory(Profile const& profile) {
    if (BoundaryConditionOption::OUTFLOW == profile.m_boundaryConditionOption) {
        return std::make_unique<OutflowBoundaryCondition>();
    }
    return nullptr;
}

} // namespace MHD
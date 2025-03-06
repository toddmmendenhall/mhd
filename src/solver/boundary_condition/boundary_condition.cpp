#include <boundary_condition/boundary_condition.hpp>
#include <execution_controller.hpp>
#include <grid.hpp>
#include <profile.hpp>
#include <variable_store.hpp>

#include <memory>

namespace MHD {

struct OutflowBoundaryConditionKernel {
    OutflowBoundaryConditionKernel(BoundaryConditionContext& context) :
        m_context(context) {}

    void operator()(std::size_t i) {
        // std::size_t const iBnd = m_context.boundaryFaceToBoundaryCellIndices[i];
        // std::size_t const iInt = m_context.boundaryFaceToInteriorCellIndices[i];

        // // Normal component of the velocity inside the boundary
        // auto uBar = m_context.u[iInt] * m_context.faceNormalX[i] +
        //             m_context.v[iInt] * m_context.faceNormalY[i] +
        //             m_context.w[iInt] * m_context.faceNormalZ[i];

        // // Neumann condition on the normal velocity
        // m_context.u[iBnd] = uBar * m_context.faceNormalX[i];
        // m_context.v[iBnd] = uBar * m_context.faceNormalY[i];
        // m_context.w[iBnd] = uBar * m_context.faceNormalZ[i];

        // // Dirichlet condition on the pressure
        // m_context.rho[iBnd] = m_context.rho[iInt];
        // m_context.p[iBnd] = m_context.p[iInt];
    }
    BoundaryConditionContext& m_context;
};

struct ReflectiveBoundaryConditionKernel {
    ReflectiveBoundaryConditionKernel(BoundaryConditionContext& context) : m_context(context) {}

    void operator()(std::size_t const i) {
        std::size_t const iIn = m_context.boundaryIdxToCellIdxs.at(m_context.boundaryIdxs[i])[0];
        std::size_t const iOut = m_context.boundaryIdxToCellIdxs.at(m_context.boundaryIdxs[i])[1];

        m_context.rho[iOut] = m_context.rho[iIn];
        m_context.p[iOut] = m_context.p[iIn];
        m_context.u[iOut] = -m_context.u[iIn];
        m_context.v[iOut] = -m_context.v[iIn];
        m_context.w[iOut] = -m_context.w[iIn];
    }

    BoundaryConditionContext& m_context;
};

BoundaryConditionContext::BoundaryConditionContext(IGrid const& grid, VariableStore& vs) :
    numBoundaries(grid.NumBoundaries()), boundaryIdxToCellIdxs(grid.BoundaryIdxToCellIdxs()),
    boundaryIdxs(grid.BoundaryIdxs()),
    faceNormalX(grid.FaceNormalX()), faceNormalY(grid.FaceNormalY()), faceNormalZ(grid.FaceNormalZ()),
    rho(vs.rho), u(vs.u), v(vs.v), w(vs.w), p(vs.p), e(vs.e), t(vs.t), cs(vs.cs) {}

class OutflowBoundaryCondition : public IBoundaryCondition {
public:
    OutflowBoundaryCondition(IGrid const& grid, VariableStore& vs) {
        m_context = std::make_unique<BoundaryConditionContext>(grid, vs);
    };

    void Compute(ExecutionController const& execCtrl) {
        OutflowBoundaryConditionKernel kern(*m_context);
        execCtrl.LaunchKernel(kern, m_context->numBoundaries);
    }
};

class ReflectiveBoundaryCondition : public IBoundaryCondition {
public:
    ReflectiveBoundaryCondition(IGrid const& grid, VariableStore& vs) {
        m_context = std::make_unique<BoundaryConditionContext>(grid, vs);
    };

    void Compute(ExecutionController const& execCtrl) {
        ReflectiveBoundaryConditionKernel kern(*m_context);
        execCtrl.LaunchKernel(kern, m_context->numBoundaries);
    }
};

std::unique_ptr<IBoundaryCondition> boundaryConditionFactory(Profile const& profile, IGrid const& grid, VariableStore& vs) {
    if (BoundaryConditionOption::REFLECTIVE == profile.m_boundaryConditionOption) {
        return std::make_unique<ReflectiveBoundaryCondition>(grid, vs);
    }
    if (BoundaryConditionOption::OUTFLOW == profile.m_boundaryConditionOption) {
        return std::make_unique<OutflowBoundaryCondition>(grid, vs);
    }
    return nullptr;
}

} // namespace MHD
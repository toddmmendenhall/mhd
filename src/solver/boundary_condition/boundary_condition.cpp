#include <boundary_condition/boundary_condition.hpp>
#include <execution_controller.hpp>
#include <grid.hpp>
#include <profile.hpp>
#include <variable_store.hpp>

#include <memory>

namespace MHD {

struct OutflowBoundaryConditionKernel {
    OutflowBoundaryConditionKernel(BoundaryConditionContext& context) : m_context(context) {}

    void operator()(std::size_t const i) {
        std::size_t const iInt = m_context.boundaryIdxToCellIdxs.at(m_context.boundaryIdxs[i])[0];
        std::size_t const iExt = m_context.boundaryIdxToCellIdxs.at(m_context.boundaryIdxs[i])[1];

        // Normal component of the velocity inside the boundary
        auto uDotN = m_context.u[iInt] * m_context.faceNormalX[i] +
                     m_context.v[iInt] * m_context.faceNormalY[i] +
                     m_context.w[iInt] * m_context.faceNormalZ[i];

        // Neumann condition on the normal velocity
        m_context.u[iExt] = uDotN * m_context.faceNormalX[i];
        m_context.v[iExt] = uDotN * m_context.faceNormalY[i];
        m_context.w[iExt] = uDotN * m_context.faceNormalZ[i];

        // Dirichlet condition on the pressure
        m_context.rho[iExt] = m_context.rho[iInt];
        m_context.p[iExt] = m_context.p[iInt];
    }
    BoundaryConditionContext& m_context;
};

struct ReflectiveBoundaryConditionKernel {
    ReflectiveBoundaryConditionKernel(BoundaryConditionContext& context) : m_context(context) {}

    void operator()(std::size_t const i) {
        std::size_t const iInt = m_context.boundaryIdxToCellIdxs.at(m_context.boundaryIdxs[i])[0];
        std::size_t const iExt = m_context.boundaryIdxToCellIdxs.at(m_context.boundaryIdxs[i])[1];

        m_context.rho[iExt] = m_context.rho[iInt];
        m_context.p[iExt] = m_context.p[iInt];
        m_context.e[iExt] = m_context.e[iInt];
        m_context.t[iExt] = m_context.t[iInt];
        m_context.cs[iExt] = m_context.cs[iInt];
        m_context.u[iExt] = -m_context.u[iInt];
        m_context.v[iExt] = -m_context.v[iInt];
        m_context.w[iExt] = -m_context.w[iInt];
        m_context.bx[iExt] = -m_context.bx[iInt];
        m_context.by[iExt] = -m_context.by[iInt];
        m_context.bz[iExt] = -m_context.bz[iInt];
    }

    BoundaryConditionContext& m_context;
};

BoundaryConditionContext::BoundaryConditionContext(IGrid const& grid, VariableStore& vs) :
    numBoundaries(grid.NumBoundaries()), boundaryIdxToCellIdxs(grid.BoundaryIdxToCellIdxs()),
    boundaryIdxs(grid.BoundaryIdxs()),
    faceNormalX(grid.FaceNormalX()), faceNormalY(grid.FaceNormalY()), faceNormalZ(grid.FaceNormalZ()),
    rho(vs.rho), u(vs.u), v(vs.v), w(vs.w), p(vs.p), e(vs.e), t(vs.t), cs(vs.cs),
    bx(vs.bx), by(vs.by), bz(vs.bz) {}

class OutflowBoundaryCondition : public IBoundaryCondition {
public:
    OutflowBoundaryCondition(IGrid const& grid, VariableStore& vs) {
        m_context = std::make_unique<BoundaryConditionContext>(grid, vs);
    };

    void ApplyBoundaryConditions(ExecutionController const& execCtrl) {
        OutflowBoundaryConditionKernel kern(*m_context);
        execCtrl.LaunchKernel(kern, m_context->numBoundaries);
    }
};

class ReflectiveBoundaryCondition : public IBoundaryCondition {
public:
    ReflectiveBoundaryCondition(IGrid const& grid, VariableStore& vs) {
        m_context = std::make_unique<BoundaryConditionContext>(grid, vs);
    };

    void ApplyBoundaryConditions(ExecutionController const& execCtrl) {
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
#pragma once

#include <grid.hpp>
#include <profile.hpp>

#include <memory>

namespace MHD {

class ExecutionController;
class VariableStore;

struct BoundaryConditionContext {
    BoundaryConditionContext(IGrid const& grid, VariableStore& vs);

    std::size_t const numBoundaries;
    std::vector<std::size_t> const boundaryIdxs;
    std::map<std::size_t, std::vector<std::size_t>> boundaryIdxToCellIdxs;

    // Properties of the face
    std::vector<double> const& faceNormalX;
    std::vector<double> const& faceNormalY;
    std::vector<double> const& faceNormalZ;

    // Primitive states
    std::vector<double>& rho;
    std::vector<double>& u;
    std::vector<double>& v;
    std::vector<double>& w;
    std::vector<double>& p;
    std::vector<double>& e;
    std::vector<double>& t;
    std::vector<double>& cs;
};

class IBoundaryCondition {
public:
    virtual ~IBoundaryCondition() = default;
    virtual void ApplyBoundaryConditions(ExecutionController const& execCtrl) = 0;
    BoundaryConditionContext const& GetContext() const { return *m_context; }

protected:
    std::unique_ptr<BoundaryConditionContext> m_context;
};

std::unique_ptr<IBoundaryCondition> boundaryConditionFactory(Profile const& profile, IGrid const& grid, VariableStore& vs);

} // namespace MHD
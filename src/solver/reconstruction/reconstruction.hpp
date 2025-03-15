#pragma once

#include <grid.hpp>
#include <memory>

namespace MHD {

class ExecutionController;
class Profile;
class VariableStore;

struct ReconstructionContext {
    ReconstructionContext(VariableStore const& vs, IGrid const& grid);

    std::size_t const numFaces;
    std::vector<std::size_t> const faceIdxs;
    std::map<std::size_t, std::vector<std::size_t>> const& faceIdxToNodeIdxs;

    // Cell-centered states
    std::vector<double> const& rho;
    std::vector<double> const& u;
    std::vector<double> const& v;
    std::vector<double> const& w;
    std::vector<double> const& p;
    std::vector<double> const& e;
    std::vector<double> const& cs;
    std::vector<double> const& bx;
    std::vector<double> const& by;
    std::vector<double> const& bz;

    // Left states
    std::vector<double> rhoLeft;
    std::vector<double> uLeft;
    std::vector<double> vLeft;
    std::vector<double> wLeft;
    std::vector<double> pLeft;
    std::vector<double> eLeft;
    std::vector<double> csLeft;
    std::vector<double> bxLeft;
    std::vector<double> byLeft;
    std::vector<double> bzLeft;

    // Right states
    std::vector<double> rhoRight;
    std::vector<double> uRight;
    std::vector<double> vRight;
    std::vector<double> wRight;
    std::vector<double> pRight;
    std::vector<double> eRight;
    std::vector<double> csRight;
    std::vector<double> bxRight;
    std::vector<double> byRight;
    std::vector<double> bzRight;
};

class IReconstruction {
public:
    virtual ~IReconstruction() = default;
    virtual void ComputeLeftRightStates(ExecutionController const& execCtrl) = 0;
    ReconstructionContext const& GetContext() const { return *m_context; }

protected:
    std::unique_ptr<ReconstructionContext> m_context;
};

std::unique_ptr<IReconstruction> reconstructionFactory(Profile const& profile, VariableStore const& varStore, IGrid const& grid);

} // namespace MHD
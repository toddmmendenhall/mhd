#pragma once

#include <execution_controller.hpp>

namespace MHD {

struct ResidualContext {
    ResidualContext(IGrid const& grid, FluxContext const& flux) :
        numCells(grid.NumCells()), cellToFaceIndices(grid.GetCellIdxToFaceIdxs()), cellSize(grid.CellSize()),
        rhoFlux(flux.rhoFlux), rhoUFlux(flux.rhoUFlux), rhoVFlux(flux.rhoVFlux),
        rhoWFlux(flux.rhoWFlux), rhoEFlux(flux.rhoEFlux), cellIdxs(grid.GetCellIdxs()) {
            rhoRes.resize(numCells, 0.0);
            rhoURes.resize(numCells, 0.0);
            rhoVRes.resize(numCells, 0.0);
            rhoWRes.resize(numCells, 0.0);
            rhoERes.resize(numCells, 0.0);
        }

    std::vector<std::size_t> const& cellIdxs;
    std::size_t const numCells;
    std::map<CellIdx, FaceIdxs> const& cellToFaceIndices;
    std::vector<double> const& cellSize;

    // Face-centered fluxes
    std::vector<double> const& rhoFlux;
    std::vector<double> const& rhoUFlux;
    std::vector<double> const& rhoVFlux;
    std::vector<double> const& rhoWFlux;
    std::vector<double> const& rhoEFlux;

    // Cell-centered residuals
    std::vector<double> rhoRes;     // mass density residual
    std::vector<double> rhoURes;    // x momentum density residual
    std::vector<double> rhoVRes;    // y momentum density residual
    std::vector<double> rhoWRes;    // z momentum density residual
    std::vector<double> rhoERes;    // total energy density residual
};

struct TransportKernel {
public:
    TransportKernel(ResidualContext& context) : m_context(context) {}

    void operator()(std::size_t const i) {
        // Get the left and right face indices for this cell
        std::size_t iLeft = m_context.cellToFaceIndices.at(i).left;
        std::size_t iRight = m_context.cellToFaceIndices.at(i).right;
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
    Residual(IGrid const& grid, FluxContext const& fc) {
        m_context = std::make_unique<ResidualContext>(grid, fc);
    }
    ~Residual() = default;
    void ComputeResidual(ExecutionController const& execCtrl) {
        TransportKernel kernel(*m_context);
        execCtrl.LaunchKernel(kernel, m_context->cellIdxs);
    }
    ResidualContext const& GetContext() const { return *m_context; }

private:
    std::unique_ptr<ResidualContext> m_context;
};

} // namespace MHD
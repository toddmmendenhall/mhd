#pragma once

#include <execution_controller.hpp>

namespace MHD {

struct ResidualContext {
    ResidualContext(IGrid const& grid, FluxContext const& flux) :
        numCells(grid.NumCells()), cellToFaceIndices(grid.CellIdxToFaceIdxs()), cellSize(grid.CellSize()),
        rhoFlux(flux.rhoFlux), rhoUFlux(flux.rhoUFlux), rhoVFlux(flux.rhoVFlux),
        rhoWFlux(flux.rhoWFlux), rhoEFlux(flux.rhoEFlux), bxFlux(flux.bxFlux), byFlux(flux.byFlux), bzFlux(flux.bzFlux) {
            rhoRes.resize(numCells, 0.0);
            rhoURes.resize(numCells, 0.0);
            rhoVRes.resize(numCells, 0.0);
            rhoWRes.resize(numCells, 0.0);
            rhoERes.resize(numCells, 0.0);
            bxRes.resize(numCells, 0.0);
            byRes.resize(numCells, 0.0);
            bzRes.resize(numCells, 0.0);
        }

    std::size_t const numCells;
    std::map<std::size_t, std::vector<std::size_t>> const& cellToFaceIndices;
    std::vector<double> const& cellSize;

    // Face-centered fluxes
    std::vector<double> const& rhoFlux;
    std::vector<double> const& rhoUFlux;
    std::vector<double> const& rhoVFlux;
    std::vector<double> const& rhoWFlux;
    std::vector<double> const& rhoEFlux;
    std::vector<double> const& bxFlux;
    std::vector<double> const& byFlux;
    std::vector<double> const& bzFlux;

    // Cell-centered residuals
    std::vector<double> rhoRes;     // mass density residual
    std::vector<double> rhoURes;    // x momentum density residual
    std::vector<double> rhoVRes;    // y momentum density residual
    std::vector<double> rhoWRes;    // z momentum density residual
    std::vector<double> rhoERes;    // total energy density residual
    std::vector<double> bxRes;      // x magnetic field residual
    std::vector<double> byRes;      // z magnetic field residual
    std::vector<double> bzRes;      // z magnetic field residual
};

struct TransportKernel {
public:
    TransportKernel(ResidualContext& context) : m_context(context) {}

    void operator()(std::size_t const i) {
        // Get the left and right face indices for this cell
        std::size_t const iLeft = m_context.cellToFaceIndices.at(i)[0];
        std::size_t const iRight = m_context.cellToFaceIndices.at(i)[1];
        double c = -1.0 / m_context.cellSize[0];

        m_context.rhoRes[i] = c * (m_context.rhoFlux[iRight] - m_context.rhoFlux[iLeft]);
        m_context.rhoURes[i] = c * (m_context.rhoUFlux[iRight] - m_context.rhoUFlux[iLeft]);
        m_context.rhoVRes[i] = c * (m_context.rhoVFlux[iRight] - m_context.rhoVFlux[iLeft]);
        m_context.rhoWRes[i] = c * (m_context.rhoWFlux[iRight] - m_context.rhoWFlux[iLeft]);
        m_context.rhoERes[i] = c * (m_context.rhoEFlux[iRight] - m_context.rhoEFlux[iLeft]);
        m_context.bxRes[i] = c * (m_context.bxFlux[iRight] - m_context.bxFlux[iLeft]);
        m_context.byRes[i] = c * (m_context.byFlux[iRight] - m_context.byFlux[iLeft]);
        m_context.bzRes[i] = c * (m_context.bzFlux[iRight] - m_context.bzFlux[iLeft]);
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
        execCtrl.LaunchKernel(kernel, m_context->numCells);
    }

    ResidualContext const& GetContext() const { return *m_context; }

private:
    std::unique_ptr<ResidualContext> m_context;
};

} // namespace MHD
#pragma once

#include <constants.hpp>
#include <grid.hpp>
#include <variable_store.hpp>

#include <vector>

namespace MHD {

struct ReconstructionContext {
    ReconstructionContext(VariableStore const& vs, IGrid const& grid) : rho(vs.rho),
        u(vs.u), v(vs.v), w(vs.w), p(vs.p), e(vs.e), faceToNodeIndices(grid.FaceToNodeIndices()),
        numFaces(grid.NumFaces()) {}

    std::size_t const numFaces;
    std::vector<std::array<std::size_t, 2>> const& faceToNodeIndices;

    // Cell-centered states
    std::vector<double> const& rho;
    std::vector<double> const& u;
    std::vector<double> const& v;
    std::vector<double> const& w;
    std::vector<double> const& p;
    std::vector<double> const& e;

    // Face-centered states
    std::vector<double> rhoFace;
    std::vector<double> uFace;
    std::vector<double> vFace;
    std::vector<double> wFace;
    std::vector<double> pFace;
    std::vector<double> eFace;
};

struct FluxContext {
    FluxContext(IGrid const& grid, ReconstructionContext const& rc) :
        numFaces(grid.NumFaces()), faceToNodeIndices(grid.FaceToNodeIndices()), faceArea(grid.FaceAreas()),
        faceNormalX(grid.FaceNormalX()), faceNormalY(grid.FaceNormalY()), faceNormalZ(grid.FaceNormalZ()),
        rhoFace(rc.rhoFace), uFace(rc.uFace), vFace(rc.vFace), wFace(rc.wFace), pFace(rc.pFace), eFace(rc.eFace) {}

    std::size_t const numFaces;
    std::vector<std::array<std::size_t, 2>> const& faceToNodeIndices;

    // properties of the faces
    std::vector<double> const& faceArea;
    std::vector<double> const& faceNormalX;
    std::vector<double> const& faceNormalY;
    std::vector<double> const& faceNormalZ;

    // Face-centered states
    std::vector<double> const& rhoFace;
    std::vector<double> const& uFace;
    std::vector<double> const& vFace;
    std::vector<double> const& wFace;
    std::vector<double> const& pFace;
    std::vector<double> const& eFace;

    // Face-centered fluxes
    std::vector<double> rhoFlux;
    std::vector<double> rhoUFlux;
    std::vector<double> rhoVFlux;
    std::vector<double> rhoWFlux;
    std::vector<double> rhoEFlux;
};

struct ResidualContext {
    ResidualContext(IGrid const& grid, FluxContext const& flux) :
        numCells(grid.NumCells()), cellToFaceIndices(grid.CellToFaceIndices()), cellSize(grid.CellSize()),
        rhoFlux(flux.rhoFlux), rhoUFlux(flux.rhoUFlux), rhoVFlux(flux.rhoVFlux),
        rhoWFlux(flux.rhoWFlux), rhoEFlux(flux.rhoEFlux) {}

    std::size_t const numCells;
    std::vector<std::array<std::size_t, 2>> const& cellToFaceIndices;
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


struct ElectricFieldContext {
    ElectricFieldContext() = default;

    // Edge variables
    std::vector<std::vector<std::size_t>> const edgeToFaceIdx;

    // Face-centered
    std::vector<double> const fluxBX;
    std::vector<double> const fluxBY;
    std::vector<double> const fluxBZ;

    // Edge-centered
    std::vector<double> eX;
    std::vector<double> eY;
    std::vector<double> eZ;
};

struct MagneticFieldContext {
    MagneticFieldContext() = default;

    double const timeStep = 0;
    std::vector<double> const cellSize;

    // Face variables
    std::vector<std::vector<std::size_t>> const faceToEdgeIdx;

    // Edge-centered
    std::vector<double> const eX;
    std::vector<double> const eY;
    std::vector<double> const eZ;

    // Face-centered
    std::vector<double> bX;
    std::vector<double> bY;
    std::vector<double> bZ;
};

struct BoundaryConditionContext {
    std::vector<std::size_t> const boundaryNodeIndices;

    // properties of the face
    std::vector<double> const faceNormalX;
    std::vector<double> const faceNormalY;
    std::vector<double> const faceNormalZ;

    // Primitives on the boundary
    std::vector<double> u;
    std::vector<double> v;
    std::vector<double> w;
    std::vector<double> p;

    // Primitives inside the boundary
    std::vector<double> const uIn;
    std::vector<double> const vIn;
    std::vector<double> const wIn;

    // Primitives outide the boundary
    double const pOut = ATMOSPHERIC_PRESSURE_STP;
};

struct IntegrationContext {
    IntegrationContext(ResidualContext const& rc, VariableStore& vs, double const tStep) : 
        tStep(tStep), numCells(rc.numCells), rhoRes(rc.rhoRes), rhoURes(rc.rhoURes), rhoVRes(rc.rhoVRes),
        rhoWRes(rc.rhoWRes), rhoERes(rc.rhoERes), rho(vs.rho), rhoU(vs.rhoU), rhoV(vs.rhoV), rhoW(vs.rhoW),
        rhoE(vs.rhoE) {}

    double const tStep;
    std::size_t const numCells;

    // Cell-centered residuals
    std::vector<double> const& rhoRes;     // mass density residual
    std::vector<double> const& rhoURes;    // x momentum density residual
    std::vector<double> const& rhoVRes;    // y momentum density residual
    std::vector<double> const& rhoWRes;    // z momentum density residual
    std::vector<double> const& rhoERes;    // total energy density residual

    // Cell-centered states
    std::vector<double>& rho;
    std::vector<double>& rhoU;
    std::vector<double>& rhoV;
    std::vector<double>& rhoW;
    std::vector<double>& rhoE;
};

} // namespace MHD
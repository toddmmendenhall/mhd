#pragma once

#include <constants.hpp>
#include <grid.hpp>
#include <variable_store.hpp>

#include <vector>

namespace MHD {

struct ReconstructionContext {
    ReconstructionContext(VariableStore const& vs, IGrid const& grid) :
        rho(vs.rho), u(vs.u), v(vs.v), w(vs.w), p(vs.p), e(vs.e), cs(vs.cs), faceToNodeIndices(grid.FaceToNodeIndices()),
        numFaces(grid.NumFaces()) {
            rhoLeft.resize(numFaces, 0.0);
            uLeft.resize(numFaces, 0.0);
            vLeft.resize(numFaces, 0.0);
            wLeft.resize(numFaces, 0.0);
            pLeft.resize(numFaces, 0.0);
            eLeft.resize(numFaces, 0.0);
            csLeft.resize(numFaces, 0.0);
            rhoRight.resize(numFaces, 0.0);
            uRight.resize(numFaces, 0.0);
            vRight.resize(numFaces, 0.0);
            wRight.resize(numFaces, 0.0);
            pRight.resize(numFaces, 0.0);
            eRight.resize(numFaces, 0.0);
            csRight.resize(numFaces, 0.0);
        }

    std::size_t const numFaces;
    std::vector<std::array<std::size_t, 4>> const& faceToNodeIndices;

    // Cell-centered states
    std::vector<double> const& rho;
    std::vector<double> const& u;
    std::vector<double> const& v;
    std::vector<double> const& w;
    std::vector<double> const& p;
    std::vector<double> const& e;
    std::vector<double> const& cs;

    // Face-centered left states
    std::vector<double> rhoLeft;
    std::vector<double> uLeft;
    std::vector<double> vLeft;
    std::vector<double> wLeft;
    std::vector<double> pLeft;
    std::vector<double> eLeft;
    std::vector<double> csLeft;

    // Face-centered right states
    std::vector<double> rhoRight;
    std::vector<double> uRight;
    std::vector<double> vRight;
    std::vector<double> wRight;
    std::vector<double> pRight;
    std::vector<double> eRight;
    std::vector<double> csRight;
};

struct FluxContext {
    FluxContext(IGrid const& grid, ReconstructionContext const& rc) :
        numFaces(grid.NumFaces()), faceToNodeIndices(grid.FaceToNodeIndices()), faceArea(grid.FaceAreas()),
        faceNormalX(grid.FaceNormalX()), faceNormalY(grid.FaceNormalY()), faceNormalZ(grid.FaceNormalZ()),
        rhoLeft(rc.rhoLeft), uLeft(rc.uLeft), vLeft(rc.vLeft), wLeft(rc.wLeft), pLeft(rc.pLeft), eLeft(rc.eLeft), csLeft(rc.csLeft),
        rhoRight(rc.rhoRight), uRight(rc.uRight), vRight(rc.vRight), wRight(rc.wRight), pRight(rc.pRight), eRight(rc.eRight), csRight(rc.csRight) {
            rhoFlux.resize(numFaces, 0.0);
            rhoUFlux.resize(numFaces, 0.0);
            rhoVFlux.resize(numFaces, 0.0);
            rhoWFlux.resize(numFaces, 0.0);
            rhoEFlux.resize(numFaces, 0.0);
        }

    std::size_t const numFaces;
    std::vector<std::array<std::size_t, 4>> const& faceToNodeIndices;

    // properties of the faces
    std::vector<double> const& faceArea;
    std::vector<double> const& faceNormalX;
    std::vector<double> const& faceNormalY;
    std::vector<double> const& faceNormalZ;

    // Face-centered left states
    std::vector<double> const& rhoLeft;
    std::vector<double> const& uLeft;
    std::vector<double> const& vLeft;
    std::vector<double> const& wLeft;
    std::vector<double> const& pLeft;
    std::vector<double> const& eLeft;
    std::vector<double> const& csLeft;

    // Face-centered right states
    std::vector<double> const& rhoRight;
    std::vector<double> const& uRight;
    std::vector<double> const& vRight;
    std::vector<double> const& wRight;
    std::vector<double> const& pRight;
    std::vector<double> const& eRight;
    std::vector<double> const& csRight;

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
        rhoWFlux(flux.rhoWFlux), rhoEFlux(flux.rhoEFlux) {
            rhoRes.resize(numCells, 0.0);
            rhoURes.resize(numCells, 0.0);
            rhoVRes.resize(numCells, 0.0);
            rhoWRes.resize(numCells, 0.0);
            rhoERes.resize(numCells, 0.0);
        }

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
    BoundaryConditionContext(IGrid const& grid, VariableStore& vs) :
        boundaryFaceToBoundaryCellIndices(grid.BoundaryFaceToBoundaryCellIndices()),
        boundaryFaceToInteriorCellIndices(grid.BoundaryFaceToInteriorCellIndices()),
        faceNormalX(grid.FaceNormalX()), faceNormalY(grid.FaceNormalY()), faceNormalZ(grid.FaceNormalZ()),
        rho(vs.rho), u(vs.u), v(vs.v), w(vs.w), p(vs.p), e(vs.e), t(vs.t), cs(vs.cs) {}

    std::vector<std::size_t> const boundaryFaceToBoundaryCellIndices;
    std::vector<std::size_t> const boundaryFaceToInteriorCellIndices;

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

struct IntegrationContext {
    IntegrationContext(ResidualContext const& rc, VariableStore& vs, double const& tStep) : 
        tStep(tStep), numCells(rc.numCells), rhoRes(rc.rhoRes), rhoURes(rc.rhoURes), rhoVRes(rc.rhoVRes),
        rhoWRes(rc.rhoWRes), rhoERes(rc.rhoERes), rho(vs.rho), rhoU(vs.rhoU), rhoV(vs.rhoV), rhoW(vs.rhoW),
        rhoE(vs.rhoE) {}

    double const& tStep;
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
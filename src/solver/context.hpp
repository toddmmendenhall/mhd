#pragma once

#include <constants.hpp>
#include <grid.hpp>
#include <variable_store.hpp>

#include <vector>

namespace MHD {

struct ReconstructionContext {
    ReconstructionContext() = default;

    std::size_t numFaces = 0;
    std::vector<std::array<std::size_t, 2>> faceToNodeIndices;

    // Left states
    std::vector<double> rhoLeft;
    std::vector<double> uLeft;
    std::vector<double> vLeft;
    std::vector<double> wLeft;
    std::vector<double> pLeft;
    std::vector<double> eLeft;

    // Right states
    std::vector<double> rhoRight;
    std::vector<double> uRight;
    std::vector<double> vRight;
    std::vector<double> wRight;
    std::vector<double> pRight;
    std::vector<double> eRight;

    // Cell-centered states
    std::vector<double> const rho;
    std::vector<double> const u;
    std::vector<double> const v;
    std::vector<double> const w;
    std::vector<double> const p;
    std::vector<double> const e;
};

struct FluxContext {
    FluxContext(IGrid const& grid, ReconstructionContext const& rc) :
        numFaces(grid.NumFaces()), faceToNodeIndices(grid.FaceToNodeIndices()), faceArea(grid.FaceAreas()),
        faceNormalX(grid.FaceNormalX()), faceNormalY(grid.FaceNormalY()), faceNormalZ(grid.FaceNormalZ()),
        rhoLeft(rc.rhoLeft), uLeft(rc.uLeft), vLeft(rc.vLeft), wLeft(rc.wLeft), pLeft(rc.pLeft), eLeft(rc.eLeft),
        rhoRight(rc.rhoRight), uRight(rc.uRight), vRight(rc.vRight), wRight(rc.wRight), pRight(rc.pRight), eRight(rc.eRight) {}

    std::size_t numFaces;
    std::vector<std::array<std::size_t, 2>> faceToNodeIndices;

    // properties of the faces
    std::vector<double> faceArea;
    std::vector<double> faceNormalX;
    std::vector<double> faceNormalY;
    std::vector<double> faceNormalZ;

    // left of the face, i.e. inside the cell
    std::vector<double> rhoLeft;
    std::vector<double> uLeft;
    std::vector<double> vLeft;
    std::vector<double> wLeft;
    std::vector<double> pLeft;
    std::vector<double> eLeft;

    // right of the face, i.e. outside the cell
    std::vector<double> rhoRight;
    std::vector<double> uRight;
    std::vector<double> vRight;
    std::vector<double> wRight;
    std::vector<double> pRight;
    std::vector<double> eRight;

    // fluxes on the face, i.e. between the cells
    std::vector<double> rhoFlux;
    std::vector<double> rhoUFlux;
    std::vector<double> rhoVFlux;
    std::vector<double> rhoWFlux;
    std::vector<double> rhoEFlux;
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
    IntegrationContext() = default;

    double const timeStep = 0.0;

    std::size_t const numVars = 0;
    std::vector<double> rhoRes;     // mass density residual
    std::vector<double> rhoURes;    // x momentum density residual
    std::vector<double> rhoVRes;    // y momentum density residual
    std::vector<double> rhoWRes;    // z momentum density residual
    std::vector<double> rhoERes;    // total energy density residual
};

} // namespace MHD
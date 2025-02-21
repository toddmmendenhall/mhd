#pragma once

#include <constants.hpp>
#include <grid.hpp>
#include <variable_store.hpp>

#include <vector>

namespace MHD {

struct FluxContext {
    FluxContext(VariableStore const& varStore, IGrid const& grid) : 
        faceArea(grid.FaceAreas()), faceNormalX(grid.FaceNormalX()), faceNormalY(grid.FaceNormalY()),
        faceNormalZ(grid.FaceNormalZ())
    {
        rhoLeft = varStore.rho;
        uLeft = varStore.u;
        vLeft = varStore.v;
        wLeft = varStore.w;
        pLeft = varStore.p;
        eLeft = varStore.e;
        bXLeft = varStore.bX;
        bYLeft = varStore.bY;
        bZLeft = varStore.bZ;
        uuLeft = varStore.uu;
        bbLeft = varStore.bb;

        rhoRight = varStore.rho;
        uRight = varStore.u;
        vRight = varStore.v;
        wRight = varStore.w;
        pRight = varStore.p;
        eRight = varStore.e;
        bXRight = varStore.bX;
        bYRight = varStore.bY;
        bZRight = varStore.bZ;
        uuRight = varStore.uu;
        bbRight = varStore.bb;
    }

    // properties of the face
    std::vector<double> const faceArea;
    std::vector<double> const faceNormalX;
    std::vector<double> const faceNormalY;
    std::vector<double> const faceNormalZ;

    // primitive variables on the left of the face, i.e. inside the cell
    std::vector<double> rhoLeft;
    std::vector<double> uLeft;
    std::vector<double> vLeft;
    std::vector<double> wLeft;
    std::vector<double> pLeft;
    std::vector<double> eLeft;
    std::vector<double> bXLeft;
    std::vector<double> bYLeft;
    std::vector<double> bZLeft;

    // primitive variables on the right of the face, i.e. outside the cell
    std::vector<double> rhoRight;
    std::vector<double> uRight;
    std::vector<double> vRight;
    std::vector<double> wRight;
    std::vector<double> pRight;
    std::vector<double> eRight;
    std::vector<double> bXRight;
    std::vector<double> bYRight;
    std::vector<double> bZRight;

    // auxiliary variables on the left of the face, i.e. inside the cell
    std::vector<double> uuLeft;
    std::vector<double> bbLeft;

    // auxiliary variables on the right of the face, i.e. outside the cell
    std::vector<double> uuRight;
    std::vector<double> bbRight;

    // fluxes on the face, i.e. between the cells
    std::vector<double> rhoFlux;
    std::vector<double> rhoUFlux;
    std::vector<double> rhoVFlux;
    std::vector<double> rhoWFlux;
    std::vector<double> rhoEFlux;
    std::vector<double> bXFlux;
    std::vector<double> bYFlux;
    std::vector<double> bZFlux;
};

struct ReconstructionContext {
    ReconstructionContext() = default;

    std::size_t numCells = 0;
    std::size_t numInteriorCells = 0;

    // Left states
    std::vector<double> rhoLeft;
    std::vector<double> rhoULeft;
    std::vector<double> rhoVLeft;
    std::vector<double> rhoWLeft;
    std::vector<double> rhoELeft;

    // Right states
    std::vector<double> rhoRight;
    std::vector<double> rhoURight;
    std::vector<double> rhoVRight;
    std::vector<double> rhoWRight;
    std::vector<double> rhoERight;

    // Cell-centered states
    std::vector<double> const rho;
    std::vector<double> const rhoU;
    std::vector<double> const rhoV;
    std::vector<double> const rhoW;
    std::vector<double> const rhoE;
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
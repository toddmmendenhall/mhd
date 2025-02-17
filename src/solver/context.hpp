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
        rhoLeft = varStore.m_rho;
        uLeft = varStore.m_u;
        vLeft = varStore.m_v;
        wLeft = varStore.m_w;
        pLeft = varStore.m_p;
        eLeft = varStore.m_e;
        bXLeft = varStore.m_bX;
        bYLeft = varStore.m_bY;
        bZLeft = varStore.m_bZ;
        uuLeft = varStore.m_uu;
        bbLeft = varStore.m_bb;

        rhoRight = varStore.m_rho;
        uRight = varStore.m_u;
        vRight = varStore.m_v;
        wRight = varStore.m_w;
        pRight = varStore.m_p;
        eRight = varStore.m_e;
        bXRight = varStore.m_bX;
        bYRight = varStore.m_bY;
        bZRight = varStore.m_bZ;
        uuRight = varStore.m_uu;
        bbRight = varStore.m_bb;
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

    std::size_t const dimension = 1;
    std::vector<double> const gridSize;
    double const timeStep = 1e-2;

    // fluxes 
    std::vector<double> const rhoFluxLeft;
    std::vector<double> const rhoUFluxLeft;
    std::vector<double> const rhoVFluxLeft;
    std::vector<double> const rhoWFluxLeft;
    std::vector<double> const rhoEFluxLeft;
    std::vector<double> const bXFluxLeft;
    std::vector<double> const bYFluxLeft;
    std::vector<double> const bZFluxLeft;

    // fluxes 
    std::vector<double> const rhoFluxRight;
    std::vector<double> const rhoUFluxRight;
    std::vector<double> const rhoVFluxRight;
    std::vector<double> const rhoWFluxRight;
    std::vector<double> const rhoEFluxRight;
    std::vector<double> const bXFluxRight;
    std::vector<double> const bYFluxRight;
    std::vector<double> const bZFluxRight;

    // cell-centered conserved variables
    std::vector<double> rho;
    std::vector<double> rhoU;
    std::vector<double> rhoV;
    std::vector<double> rhoW;
    std::vector<double> rhoE;
    std::vector<double> bX;
    std::vector<double> bY;
    std::vector<double> bZ;
};

struct ElectricFieldContext {
    ElectricFieldContext() = default;

    // Edge variables
    std::vector<std::vector<std::size_t>> const m_edgeToFaceIdx;

    // Face-centered
    std::vector<double> const m_fluxBX;
    std::vector<double> const m_fluxBY;
    std::vector<double> const m_fluxBZ;

    // Edge-centered
    std::vector<double> m_eX;
    std::vector<double> m_eY;
    std::vector<double> m_eZ;
};

struct MagneticFieldContext {
    MagneticFieldContext() = default;

    double const m_timeStep = 0;
    std::vector<double> const m_cellSize;

    // Face variables
    std::vector<std::vector<std::size_t>> const m_faceToEdgeIdx;

    // Edge-centered
    std::vector<double> const m_eX;
    std::vector<double> const m_eY;
    std::vector<double> const m_eZ;

    // Face-centered
    std::vector<double> m_bX;
    std::vector<double> m_bY;
    std::vector<double> m_bZ;
};

struct BoundaryConditionContext {
    std::vector<std::size_t> const m_boundaryNodeIndices;

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

} // namespace MHD
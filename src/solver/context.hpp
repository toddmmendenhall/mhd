#pragma once

#include <vector>

namespace MHD {

struct FluxContext {
    FluxContext();

    // properties of the face
    std::vector<double> const& faceArea;
    std::vector<double> const& faceNormalX;
    std::vector<double> const& faceNormalY;
    std::vector<double> const& faceNormalZ;

    // primitive variables on the left of the face, i.e. inside the cell
    std::vector<double> const& rhoLeft;
    std::vector<double> const& uLeft;
    std::vector<double> const& vLeft;
    std::vector<double> const& wLeft;
    std::vector<double> const& pLeft;
    std::vector<double> const& eLeft;
    std::vector<double> const& bXLeft;
    std::vector<double> const& bYLeft;
    std::vector<double> const& bZLeft;

    // primitive variables on the right of the face, i.e. outside the cell
    std::vector<double> const& rhoRight;
    std::vector<double> const& uRight;
    std::vector<double> const& vRight;
    std::vector<double> const& wRight;
    std::vector<double> const& pRight;
    std::vector<double> const& eRight;
    std::vector<double> const& bXRight;
    std::vector<double> const& bYRight;
    std::vector<double> const& bZRight;

    // auxiliary variables on the left of the face, i.e. inside the cell
    std::vector<double> const& uuLeft;
    std::vector<double> const& bbLeft;

    // auxiliary variables on the right of the face, i.e. outside the cell
    std::vector<double> const& uuRight;
    std::vector<double> const& bbRight;

    // fluxes on the face, i.e. between the cells
    std::vector<double>& rhoFlux;
    std::vector<double>& rhoUFlux;
    std::vector<double>& rhoVFlux;
    std::vector<double>& rhoWFlux;
    std::vector<double>& rhoEFlux;
    std::vector<double>& bXFlux;
    std::vector<double>& bYFlux;
    std::vector<double>& bZFlux;
};

struct ReconstructionContext {
    ReconstructionContext();

    std::size_t const dimension;
    std::vector<double> const& gridSize;
    double const timeStep;

    // fluxes 
    std::vector<double> const& rhoFluxLeft;
    std::vector<double> const& rhoUFluxLeft;
    std::vector<double> const& rhoVFluxLeft;
    std::vector<double> const& rhoWFluxLeft;
    std::vector<double> const& rhoEFluxLeft;
    std::vector<double> const& bXFluxLeft;
    std::vector<double> const& bYFluxLeft;
    std::vector<double> const& bZFluxLeft;

    // fluxes 
    std::vector<double> const& rhoFluxRight;
    std::vector<double> const& rhoUFluxRight;
    std::vector<double> const& rhoVFluxRight;
    std::vector<double> const& rhoWFluxRight;
    std::vector<double> const& rhoEFluxRight;
    std::vector<double> const& bXFluxRight;
    std::vector<double> const& bYFluxRight;
    std::vector<double> const& bZFluxRight;

    // cell-centered conserved variables
    std::vector<double>& rho;
    std::vector<double>& rhoU;
    std::vector<double>& rhoV;
    std::vector<double>& rhoW;
    std::vector<double>& rhoE;
    std::vector<double>& bX;
    std::vector<double>& bY;
    std::vector<double>& bZ;
};

struct ElectricFieldContext {
    ElectricFieldContext();

    // Edge variables
    std::vector<std::vector<std::size_t>> const& m_edgeToFaceIdx;

    // Face-centered
    std::vector<double> const& m_fluxBX;
    std::vector<double> const& m_fluxBY;
    std::vector<double> const& m_fluxBZ;

    // Edge-centered
    std::vector<double>& m_eX;
    std::vector<double>& m_eY;
    std::vector<double>& m_eZ;
};

struct MagneticFieldContext {
    MagneticFieldContext();

    double const m_timeStep;
    std::vector<double> const& m_cellSize;

    // Face variables
    std::vector<std::vector<std::size_t>> const& m_faceToEdgeIdx;

    // Edge-centered
    std::vector<double> const& m_eX;
    std::vector<double> const& m_eY;
    std::vector<double> const& m_eZ;

    // Face-centered
    std::vector<double>& m_bX;
    std::vector<double>& m_bY;
    std::vector<double>& m_bZ;
};

} // namespace MHD
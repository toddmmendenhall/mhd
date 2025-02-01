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

} // namespace MHD
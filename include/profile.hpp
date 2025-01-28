#pragma once

#include <profile_options.hpp>

#include <vector>

namespace MHD {

class Profile {
public:
    Profile();
    ~Profile() = default;

    Dimension const GetGridDimension() const;
    Geometry const GetGridGeometry() const;
    std::vector<double> const& GetGridBounds() const;
    std::vector<double> const& GetGridSpacing() const;
    BoundaryCondition const& GetGridBoundaryConditions() const;
    SpatialDerivativeMethod const GetSpatialDerivativeMethod() const;
    TemporalIntegrationMethod const GetTemporalIntegrationMethod() const;
    FluxScheme const GetFluxScheme() const {return m_fluxScheme;}

    void SetGridDimension(Dimension const dimension);
    void SetGridGeometry(Geometry const geometry);
    void SetGridBounds(std::vector<double> const& bounds);
    void SetGridSpacing(std::vector<double> const& spacing);
    void SetGridBoundaryConditions(BoundaryCondition const& boundaryConditions);
    void SetSpatialDerivativeMethod(SpatialDerivativeMethod const method);
    void SetTemporalIntegrationMethod(TemporalIntegrationMethod const method);
    void SetFluxScheme(FluxScheme const value) {m_fluxScheme = value;}

private:
    Dimension m_gridDimension;
    Geometry m_gridGeometry;
    std::vector<double> m_gridBounds;
    std::vector<double> m_gridSpacing;
    BoundaryCondition m_gridBoundaryConditions;
    SpatialDerivativeMethod m_spatialDerivativeMethod;
    TemporalIntegrationMethod m_temporalIntegrationMethod;
    FluxScheme m_fluxScheme;
};

} // namespace MHD

#pragma once

#include <variable_store.hpp>

#include <cmath>
#include <cstddef>
#include <vector>

namespace MHD {

struct VelocityKernel {
    VelocityKernel(VariableStore& vs) :
        rho(vs.rho), rhoU(vs.rhoU), rhoV(vs.rhoV), rhoW(vs.rhoW), u(vs.u), v(vs.v), w(vs.w) {}
    
    inline void operator()(std::size_t const i) {
        if (rho[i] < 0.0) {
            throw;
        }
        double const rhoInv = 1.0 / rho[i];
        u[i] = rhoU[i] * rhoInv;
        v[i] = rhoV[i] * rhoInv;
        w[i] = rhoW[i] * rhoInv;
    }

    std::vector<double> const& rho;
    std::vector<double> const& rhoU;
    std::vector<double> const& rhoV;
    std::vector<double> const& rhoW;
    std::vector<double>& u;
    std::vector<double>& v;
    std::vector<double>& w;
};

struct MagneticFieldSquaredKernel {
    MagneticFieldSquaredKernel(std::vector<double> const& bx, std::vector<double> const& by, std::vector<double> const& bz,
                               std::vector<double>& b_b) :
        bx(bx), by(by), bz(bz), b_b(b_b) {}
    
    inline void operator()(std::size_t idx) {
        b_b[idx] = bx[idx] * bx[idx] + by[idx] * by[idx] + bz[idx] * bz[idx];
    }

    std::vector<double> const& bx;
    std::vector<double> const& by;
    std::vector<double> const& bz;
    std::vector<double>& b_b;
};

struct SpecificInternalEnergyKernel {
    SpecificInternalEnergyKernel(VariableStore& vs) :
        rho(vs.rho), rhoE(vs.rhoE), rhoU(vs.rhoU), rhoV(vs.rhoV), rhoW(vs.rhoW), e(vs.e) {}

    inline void operator()(std::size_t const i) {
        double const rhoInv = 1.0 / rho[i];
        e[i] = (rhoE[i] - 0.5 * (rhoU[i] * rhoU[i] + rhoV[i] * rhoV[i] + rhoW[i] * rhoW[i]) * rhoInv) * rhoInv;
    }

    std::vector<double> const& rho;
    std::vector<double> const& rhoE;
    std::vector<double> const& rhoU;
    std::vector<double> const& rhoV;
    std::vector<double> const& rhoW;
    std::vector<double>& e;
};

struct CaloricallyPerfectGasPressureKernel {
    CaloricallyPerfectGasPressureKernel(VariableStore& vs) :
        gammaMinusOne(vs.gamma - 1.0), rho(vs.rho), e(vs.e), p(vs.p) {}

    inline void operator()(std::size_t const i) {
        if (rho[i] < 0.0 || e[i] < 0.0) {
            throw;
        }
        p[i] = gammaMinusOne * rho[i] * e[i];
    }

    double const gammaMinusOne;
    std::vector<double> const& rho;
    std::vector<double> const& e;
    std::vector<double>& p;
};

struct CaloricallyPerfectGasTemperatureKernel {
    CaloricallyPerfectGasTemperatureKernel(VariableStore& vs) :
        gammaMinusOne(vs.gamma - 1.0), rInv(1.0 / vs.r), e(vs.e), t(vs.t) {}

    inline void operator()(std::size_t const i) {
        t[i] = gammaMinusOne * rInv * e[i];
    }

    double const gammaMinusOne;
    double const rInv;
    std::vector<double> const& e;
    std::vector<double>& t;
};

struct CaloricallyPerfectGasSoundSpeedKernel {
    CaloricallyPerfectGasSoundSpeedKernel(VariableStore& vs) :
        gammaTimesGammaMinusOne(vs.gamma * (vs.gamma - 1.0)), e(vs.e), cs(vs.cs) {}

    inline void operator()(std::size_t const i) {
        cs[i] = std::sqrt(gammaTimesGammaMinusOne * e[i]);
    }

    double const gammaTimesGammaMinusOne;
    std::vector<double> const& e;
    std::vector<double>& cs;
};

struct MaximumWaveSpeedKernel {
    MaximumWaveSpeedKernel(VariableStore& vs) :
        u(vs.u), cs(vs.cs), sMax(vs.sMax) { sMax = 0.0; }

    inline void operator()(std::size_t const i) {
        double const waveSpeed = std::abs(u[i]) + cs[i];
        if (waveSpeed > sMax) {
            sMax = waveSpeed;
        }
    }

    std::vector<double> const& u;
    std::vector<double> const& cs;
    double& sMax;
};

struct MomentumDensityKernel {
    MomentumDensityKernel(VariableStore& vs) :
        rho(vs.rho), u(vs.u), v(vs.v), w(vs.w), rhoU(vs.rhoU), rhoV(vs.rhoV), rhoW(vs.rhoW) {}
    
    void operator()(std::size_t const i) {
        rhoU[i] = rho[i] * u[i];
        rhoV[i] = rho[i] * v[i];
        rhoW[i] = rho[i] * w[i];
    }

    std::vector<double> const& rho;
    std::vector<double> const& u;
    std::vector<double> const& v;
    std::vector<double> const& w;
    std::vector<double>& rhoU;
    std::vector<double>& rhoV;
    std::vector<double>& rhoW;
};

struct TotalEnergyDensityKernel {
    TotalEnergyDensityKernel(VariableStore& vs) :
        rho(vs.rho), u(vs.u), v(vs.v), w(vs.w), e(vs.e), rhoE(vs.rhoE) {}

    void operator()(std::size_t const i) {
        rhoE[i] = rho[i] * (e[i] + 0.5 * (u[i] * u[i] + v[i] * v[i] + w[i] * w[i]));
    }

    std::vector<double> const& rho;
    std::vector<double> const& u;
    std::vector<double> const& v;
    std::vector<double> const& w;
    std::vector<double> const& e;
    std::vector<double>& rhoE;
};

} // namespace MHD
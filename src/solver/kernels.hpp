#pragma once

#include <variable_store.hpp>

#include <cmath>
#include <cstddef>
#include <vector>

namespace MHD {

struct SpecificVolumeKernel {
    SpecificVolumeKernel(VariableStore& vs) :
        rho(vs.rho), rhoInv(vs.rhoInv) {}

    inline void operator()(std::size_t idx) {
        rhoInv[idx] = 1.0 / rho[idx];
    }

    std::vector<double> const& rho;
    std::vector<double>& rhoInv;
};

struct VelocityKernel {
    VelocityKernel(VariableStore& vs) :
        rhoInv(vs.rhoInv), rhoU(vs.rhoU), rhoV(vs.rhoV), rhoW(vs.rhoW), u(vs.u), v(vs.v), w(vs.w), uu(vs.uu) {}
    
    inline void operator()(std::size_t idx) {
        u[idx] = rhoU[idx] * rhoInv[idx];
        v[idx] = rhoV[idx] * rhoInv[idx];
        w[idx] = rhoW[idx] * rhoInv[idx];
        uu[idx] = u[idx] * u[idx] + v[idx] * v[idx] + w[idx] * w[idx];
    }

    std::vector<double> const& rhoInv;
    std::vector<double> const& rhoU;
    std::vector<double> const& rhoV;
    std::vector<double> const& rhoW;
    std::vector<double>& u;
    std::vector<double>& v;
    std::vector<double>& w;
    std::vector<double>& uu;
};

struct VelocitySquaredKernel {
    VelocitySquaredKernel(VariableStore& vs) :
        u(vs.u), v(vs.v), w(vs.w), uu(vs.uu) {}
    
    inline void operator()(std::size_t idx) {
        uu[idx] = u[idx] * u[idx] + v[idx] * v[idx] + w[idx] * w[idx];
    }

    std::vector<double> const& u;
    std::vector<double> const& v;
    std::vector<double> const& w;
    std::vector<double>& uu;
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
        rhoInv(vs.rhoInv), rhoE(vs.rhoE), uu(vs.uu), e(vs.e) {}

    inline void operator()(std::size_t idx) {
        e[idx] = rhoE[idx] * rhoInv[idx] - 0.5 * uu[idx];
        if (e[idx] < 0.0) {
            e[idx] = 0.0;
        }
    }

    std::vector<double> const& rhoInv;
    std::vector<double> const& rhoE;
    std::vector<double> const& uu;
    std::vector<double>& e;
};

struct CaloricallyPerfectGasPressureKernel {
    CaloricallyPerfectGasPressureKernel(VariableStore& vs) :
        gammaMinus1(vs.gamma - 1.0), rho(vs.rho), e(vs.e), p(vs.p) {}

    inline void operator()(std::size_t idx) {
        if (rho[idx] < 1e-9 || e[idx] < 1e-9) {
            p[idx] = 0.0;
        }
        p[idx] = gammaMinus1 * rho[idx] * e[idx];
    }

    double const gammaMinus1;
    std::vector<double> const& rho;
    std::vector<double> const& e;
    std::vector<double>& p;
};

struct CaloricallyPerfectGasSpecificInternalEnergyKernel {
    CaloricallyPerfectGasSpecificInternalEnergyKernel(VariableStore& vs) :
        gammaMinus1Inv(1.0 / (vs.gamma - 1.0)), rhoInv(vs.rhoInv), p(vs.p), e(vs.e) {}

    inline void operator()(std::size_t i) {
        e[i] = p[i] * gammaMinus1Inv * rhoInv[i];
    }

    double const gammaMinus1Inv;
    std::vector<double> const& rhoInv;
    std::vector<double> const& p;
    std::vector<double>& e;
};

struct PerfectGasTemperatureKernel {
    PerfectGasTemperatureKernel(VariableStore& vs) :
        rInv(1.0 / vs.r), rhoInv(vs.rhoInv), p(vs.p), t(vs.t) {}

    inline void operator()(std::size_t idx) {
        t[idx] = p[idx] * rhoInv[idx] * rInv;
    }

    double const rInv;
    std::vector<double> const& rhoInv;
    std::vector<double> const& p;
    std::vector<double>& t;
};

double inline CpOverR(std::vector<double> const& a, double const T) {
    double const T2 = T * T;
    return a[0] / T2 + a[1] / T + a[2] + a[3] * T + a[4] * T2 + a[5] * T2 * T + a[6] * T2 * T2;
}

double inline HOverRT(std::vector<double> const& a, double const T) {
    double const T2 = T * T;
    return -a[0] / T2 + a[1] * std::log(T) / T + a[2] + a[3] / 2 * T + a[4] / 3 * T2 + a[5] / 4 * T2 * T + a[6] / 5 * T2 * T2 + a[7] / T;
}

double inline SOverR(std::vector<double> const& a, double const T) {
    double const T2 = T * T;
    return -a[0] / 2 / T2 - a[1] / T + a[2] * std::log(T) + a[3] * T + a[4]/ 2 * T2 + a[5] / 3 * T2 * T + a[6] / 4 * T2 * T2 + a[8];
}

struct ThermallyPerfectGasPressureKernel {
    ThermallyPerfectGasPressureKernel(std::vector<double> const& rho, std::vector<double>& p) :
        m_rho(rho), m_p(p) {}

    void operator()(std::size_t idx) {
        m_p[idx] = m_rho[idx];
    }

    std::vector<double> const& m_rho;
    std::vector<double>& m_p;
};

struct MomentumDensityKernel {
    MomentumDensityKernel(VariableStore& vs) :
        rho(vs.rho), u(vs.u), v(vs.v), w(vs.w), rhoU(vs.rhoU), rhoV(vs.rhoV), rhoW(vs.rhoW) {}
    
    void operator()(std::size_t i) {
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
        rho(vs.rho), uu(vs.uu), e(vs.e), rhoE(vs.rhoE) {}

    void operator()(std::size_t idx) {
        rhoE[idx] = rho[idx] * e[idx] + 0.5 * uu[idx];
    }

    std::vector<double> const& rho;
    std::vector<double> const& uu;
    std::vector<double> const& e;
    std::vector<double>& rhoE;
};

} // namespace MHD
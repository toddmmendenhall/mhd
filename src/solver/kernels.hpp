#pragma once

#include <constants.hpp>

#include <cmath>
#include <cstddef>
#include <vector>

namespace MHD {

struct SpecificVolumeKernel {
    SpecificVolumeKernel(std::vector<double> const& rho,
                         std::vector<double>& rhoInv) :
        rho(rho), rhoInv(rhoInv) {}

    inline void operator()(std::size_t idx) {
        rhoInv[idx] = 1.0 / rho[idx];
    }

    std::vector<double> const& rho;
    std::vector<double>& rhoInv;
};

struct VelocityKernel {
    VelocityKernel(std::vector<double> const& rhoInv,
                   std::vector<double> const& rho_u, std::vector<double> const& rho_v, std::vector<double> const& rho_w,
                   std::vector<double>& u, std::vector<double>& v, std::vector<double>& w, std::vector<double>& u_u) :
        rhoInv(rhoInv), rho_u(rho_u), rho_v(rho_v), rho_w(rho_w), u(u), v(v), w(w), u_u(u_u) {}
    
    inline void operator()(std::size_t idx) {
        u[idx] = rho_u[idx] * rhoInv[idx];
        v[idx] = rho_v[idx] * rhoInv[idx];
        w[idx] = rho_w[idx] * rhoInv[idx];
        u_u[idx] = u[idx] * u[idx] + v[idx] * v[idx] + w[idx] * w[idx];
    }

    std::vector<double> const& rhoInv;
    std::vector<double> const& rho_u;
    std::vector<double> const& rho_v;
    std::vector<double> const& rho_w;
    std::vector<double>& u;
    std::vector<double>& v;
    std::vector<double>& w;
    std::vector<double>& u_u;
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
    SpecificInternalEnergyKernel(std::vector<double> const& rhoInv, std::vector<double> const& rhoE, std::vector<double> const& uu,
                                 std::vector<double> const& bb, std::vector<double>& e) :
        m_rhoInv(rhoInv), m_rhoE(rhoE), m_uu(uu), m_bb(bb), m_e(e) {}

    inline void operator()(std::size_t idx) {
        m_e[idx] = m_rhoE[idx] * m_rhoInv[idx] - 0.5 * m_uu[idx] - 0.5 * VACUUM_PERMEABILITY_INV * m_bb[idx] * m_rhoInv[idx];
    }

    std::vector<double> const& m_rhoInv;
    std::vector<double> const& m_rhoE;
    std::vector<double> const& m_uu;
    std::vector<double> const& m_bb;
    std::vector<double>& m_e;
};

struct CaloricallyPerfectGasPressureKernel {
    CaloricallyPerfectGasPressureKernel(double const gamma,
                                        std::vector<double> const& rho, std::vector<double> const& e, std::vector<double> const& bb,
                                        std::vector<double>& p) :
        m_gammaMinus1(gamma - 1.0), m_rho(rho), m_e(e), m_bb(bb), m_p(p) {}

    inline void operator()(std::size_t idx) {
        m_p[idx] = m_rho[idx] * m_gammaMinus1 * m_e[idx] + 0.5 * VACUUM_PERMEABILITY_INV * m_bb[idx];
    }

    double const m_gammaMinus1;
    std::vector<double> const& m_rho;
    std::vector<double> const& m_e;
    std::vector<double> const& m_bb;
    std::vector<double>& m_p;
};

struct PerfectGasTemperatureKernel {
    PerfectGasTemperatureKernel(double const r, std::vector<double> const& rhoInv,
                                std::vector<double> const& p, std::vector<double>& t) :
        m_rInv(1.0 / r), m_rhoInv(rhoInv), m_p(p), m_t(t) {}

    inline void operator()(std::size_t idx) {
        m_t[idx] = m_p[idx] * m_rhoInv[idx] * m_rInv;
    }

    double const m_rInv;
    std::vector<double> const& m_rhoInv;
    std::vector<double> const& m_p;
    std::vector<double>& m_t;
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
    MomentumDensityKernel(std::vector<double> const& rho,
                          std::vector<double> const& u, std::vector<double> const& v, std::vector<double> const& w,
                          std::vector<double>& rhoU, std::vector<double>& rhoV, std::vector<double>& rhoW) :
        m_rho(rho), m_u(u), m_v(v), m_w(w), m_rhoU(rhoU), m_rhoV(rhoV), m_rhoW(rhoW) {}
    
    void operator()(std::size_t idx) {
        m_rhoU[idx] = m_rho[idx] * m_u[idx];
        m_rhoV[idx] = m_rho[idx] * m_v[idx];
        m_rhoW[idx] = m_rho[idx] * m_w[idx];
    }

    std::vector<double> const& m_rho;
    std::vector<double> const& m_u;
    std::vector<double> const& m_v;
    std::vector<double> const& m_w;
    std::vector<double>& m_rhoU;
    std::vector<double>& m_rhoV;
    std::vector<double>& m_rhoW;
};

struct TotalEnergyDensityKernel {
    TotalEnergyDensityKernel(std::vector<double> const& rho, std::vector<double> const& rhoInv,
                             std::vector<double> const& uu, std::vector<double> const& e,
                             std::vector<double> const& bb, std::vector<double>& rhoE) :
        m_rho(rho), m_rhoInv(rhoInv), m_uu(uu), m_e(e), m_bb(bb), m_rhoE(rhoE) {}

    void operator()(std::size_t idx) {
        m_rhoE[idx] = m_rho[idx] * m_e[idx] + 0.5 * m_uu[idx] + 0.5 * VACUUM_PERMEABILITY_INV * m_bb[idx] * m_rhoInv[idx];
    }

    std::vector<double> const& m_rho;
    std::vector<double> const& m_rhoInv;
    std::vector<double> const& m_uu;
    std::vector<double> const& m_e;
    std::vector<double> const& m_bb;
    std::vector<double>& m_rhoE;
};

} // namespace MHD
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
    SpecificInternalEnergyKernel(std::vector<double> const& rhoInv, std::vector<double> const& rho_e, std::vector<double> const& u_u,
                                 std::vector<double> const& b_b, std::vector<double>& eInt) :
        rhoInv(rhoInv), rho_e(rho_e), u_u(u_u), b_b(b_b), eInt(eInt) {}

    inline void operator()(std::size_t idx) {
        eInt[idx] = rho_e[idx] * rhoInv[idx] - 0.5 * u_u[idx] - 0.5 * VACUUM_PERMITTIVITY_INV * b_b[idx] * rhoInv[idx];
    }

    std::vector<double> const& rhoInv;
    std::vector<double> const& rho_e;
    std::vector<double> const& u_u;
    std::vector<double> const& b_b;
    std::vector<double>& eInt;
};

struct CaloricallyPerfectGasPressureKernel {
    CaloricallyPerfectGasPressureKernel(double const gammaMinus1,
                                        std::vector<double> const& rho, std::vector<double> const& eInt, std::vector<double> const& b_b,
                                        std::vector<double>& pres) :
        gammaMinus1(gammaMinus1), rho(rho), eInt(eInt), b_b(b_b), pres(pres) {}

    inline void operator()(std::size_t idx) {
        pres[idx] = rho[idx] * gammaMinus1 * eInt[idx] + 0.5 * VACUUM_PERMITTIVITY_INV * b_b[idx];
    }

    double const gammaMinus1;
    std::vector<double> const& rho;
    std::vector<double> const& eInt;
    std::vector<double> const& b_b;
    std::vector<double>& pres;
};

struct PerfectGasTemperatureKernel {
    PerfectGasTemperatureKernel(double const rSpecInv, std::vector<double> const& rhoInv,
                                std::vector<double> const& pres, std::vector<double>& temp) :
        rSpecInv(rSpecInv), rhoInv(rhoInv), pres(pres), temp(temp) {}

    inline void operator()(std::size_t idx) {
        temp[idx] = pres[idx] * rhoInv[idx] * rSpecInv;
    }

    double const rSpecInv;
    std::vector<double> const& rhoInv;
    std::vector<double> const& pres;
    std::vector<double>& temp;
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
        rho(rho), p(p) {}

    void operator()(std::size_t idx) {
        p[idx] = rho[idx];
    }

    std::vector<double> const& rho;
    std::vector<double>& p;
};

} // namespace MHD
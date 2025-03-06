#pragma once

#include <variable_store.hpp>

#pragma once

#include <variable_store.hpp>

#include <cmath>
#include <cstddef>
#include <vector>

namespace MHD {

struct CaloricallyPerfectGasPressureKernel {
    CaloricallyPerfectGasPressureKernel(VariableStore& vs) :
        gammaMinusOne(vs.gamma - 1.0), rho(vs.rho), e(vs.e), p(vs.p) {}

    inline void operator()(std::size_t const i) {
        if (rho[i] < 1e-9 || e[i] < 1e-9) {
            p[i] = 0.0;
        }
        p[i] = gammaMinusOne * rho[i] * e[i];
    }

    double const gammaMinusOne;
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

struct CaloricallyPerfectGasSoundSpeedKernel {
    CaloricallyPerfectGasSoundSpeedKernel(VariableStore& vs) :
        r(vs.r), gamma(vs.gamma), t(vs.t), cs(vs.cs) {}

    inline void operator()(std::size_t const i) {
        cs[i] = std::sqrt(gamma * r * t[i]);
    }

    double const r;
    double const gamma;
    std::vector<double> const& t;
    std::vector<double>& cs;
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

} // namespace MHD
#include <cmath>
#include <cstddef>
#include <vector>

struct VelocityUKernel {
    VelocityUKernel(std::vector<double> const& rho, std::vector<double> const& rho_u, std::vector<double>& u) :
        rho(rho), rho_u(rho_u), u(u) {}
    
    void operator()(std::size_t idx) {
        u[idx] = rho_u[idx] / rho[idx];
    }

    std::vector<double> const& rho;
    std::vector<double> const& rho_u;
    std::vector<double>& u;
};

struct InternalEnergyKernel {
    InternalEnergyKernel(std::vector<double> const& rho, std::vector<double> const& rho_u, std::vector<double> const& rho_e0, std::vector<double>& e) :
        rho(rho), rho_u(rho_u), rho_e0(rho_e0), e(e) {}

    void operator()(std::size_t idx) {
        double const rhoInverse = 1.0 / rho[idx];
        double const u = rho_u[idx] * rhoInverse;
        e[idx] = rho_e0[idx] * rhoInverse - 0.5 * u * u;
    }

    std::vector<double> const& rho;
    std::vector<double> const& rho_u;
    std::vector<double> const& rho_e0;
    std::vector<double>& e;
};

struct CaloricallyPerfectIdealGasPressureKernel {
    CaloricallyPerfectIdealGasPressureKernel(std::vector<double> const& rho, std::vector<double> const& rho_u, std::vector<double> const& rho_e0,
                                             double const gammaMinusOne, std::vector<double>& P) :
        rho(rho), rho_u(rho_u), rho_e0(rho_e0), gammaMinusOne(gammaMinusOne), P(P) {}

    void operator()(std::size_t idx) {
        P[idx] = gammaMinusOne * (rho_e0[idx] - 0.5 * rho_u[idx] * rho_u[idx] / rho[idx]);
    }

    std::vector<double> const& rho;
    std::vector<double> const& rho_u;
    std::vector<double> const& rho_e0;
    double const gammaMinusOne;
    std::vector<double>& P;
};

struct PerfectGasTemperatureKernel {
    PerfectGasTemperatureKernel(std::vector<double> const& rho, std::vector<double> const& pressure, double const rSpecific, std::vector<double>& T) :
        rho(rho), pressure(pressure), R(rSpecific), T(T) {}

    void operator()(std::size_t idx) {
        T[idx] = pressure[idx] / (rho[idx] * R);
    }

    std::vector<double> const& rho;
    std::vector<double> const& pressure;
    double const R;
    std::vector<double>& T;
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

    double operator()(std::size_t idx) {
        p[idx] = rho[idx];
    }

    std::vector<double> const& rho;
    std::vector<double>& p;
};
#pragma once

namespace MHD {

double constexpr ATMOSPHERIC_DENSITY_STP = 1.204;                       // kg / m^3
double constexpr ATMOSPHERIC_PRESSURE_STP = 101325;                     // Pa

double constexpr VACUUM_PERMEABILITY = 1.256E-6;                        // N * A^-2
double constexpr VACUUM_PERMEABILITY_INV = 1.0 / VACUUM_PERMEABILITY;   // A^2 / N

double constexpr PI = 3.141592654;
double constexpr TWO_PI = 2.0 * PI;
double constexpr FOUR_PI = 4.0 * PI;

} // namespace MHD

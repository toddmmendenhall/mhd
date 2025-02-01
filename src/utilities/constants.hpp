#pragma once

namespace MHD {

double constexpr ATMOSPHERIC_DENSITY_STP = 1.204;                       // kg / m^3
double constexpr ATMOSPHERIC_PRESSURE_STP = 101325;                     // Pa

double constexpr VACUUM_PERMITTIVITY = 1.256E-6;                        // N * A^-2
double constexpr VACUUM_PERMITTIVITY_INV = 1.0 / VACUUM_PERMITTIVITY;   // A^2 / N

} // namespace MHD

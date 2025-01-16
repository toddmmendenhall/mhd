#include <map>
#include <string>
#include <vector>

struct SpeciesData {
    SpeciesData() = default;
    SpeciesData(double const molecularMass,
                std::size_t const numTemperatureIntervals,
                std::size_t const numCoefficientsPerInterval,
                std::vector<double> const& temperatures,
                std::vector<double> const& coefficients) :
        m_molecularMass(molecularMass),
        m_numTemperatureIntervals(numTemperatureIntervals),
        m_numCoefficientsPerInterval(numCoefficientsPerInterval),
        m_temperatures(temperatures),
        m_coefficients(coefficients) {}

    double m_molecularMass; // g/mol
    std::size_t m_numTemperatureIntervals;
    std::size_t m_numCoefficientsPerInterval;
    std::vector<double> m_temperatures;  // K
    std::vector<double> m_coefficients;
};

struct ThermodynamicsData {
    ThermodynamicsData();
    std::map<std::string, SpeciesData> m_speciesData;
};
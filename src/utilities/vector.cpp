#include "vector.hpp"

#include <cmath>

namespace GEOM_UTILS {

Vector2D::Vector2D(double const x, double const y) {
    m_components[0] = x;
    m_components[1] = y;
}

Vector3D::Vector3D(double x, double y, double z) {
    m_components[0] = x;
    m_components[1] = y;
    m_components[2] = z;
}

double const Vector3D::X() const {
    return m_components[0];
}

double const Vector3D::Y() const {
    return m_components[1];
}

double const Vector3D::Z() const {
    return m_components[2];
}

void Vector3D::X(double const c) {
    m_components[0] = c;
}

void Vector3D::Y(double const c) {
    m_components[1] = c;
}

void Vector3D::Z(double const c) {
    m_components[2] = c;
}

double const Vector3D::MagnitudeSquared() const {
    return m_components[0] * m_components[0] + m_components[1] * m_components[1] + m_components[2] * m_components[2];
}

double const Vector3D::Magnitude() const {
    return std::sqrt(MagnitudeSquared());
}

} // namespace GEOM_UTILS

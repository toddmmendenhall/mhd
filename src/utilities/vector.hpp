#pragma once

#include <array>
#include <cmath>

class Vec3 {
private:
    std::array<double, 3> m_components = {NAN, NAN, NAN};

public:
    Vec3() = default;

    Vec3(double const x, double const y, double const z) {
        m_components[0] = x;
        m_components[1] = y;
        m_components[2] = z;
    }

    Vec3(std::array<double, 3> const& v) {
        m_components[0] = v[0];
        m_components[1] = v[1];
        m_components[2] = v[2];
    }

    ~Vec3() = default;

    double Norm() const {
        return std::sqrt(NormSquared());
    }

    double NormSquared() const {
        return m_components[0] * m_components[0] + m_components[1] * m_components[1] + m_components[2] * m_components[2];
    }

    double operator[](std::size_t const i) const {
        return m_components[i];
    }

    double X() const {
        return m_components[0];
    }

    double Y() const {
        return m_components[1];
    }

    double Z() const {
        return m_components[2];
    }

    double& operator[](std::size_t const i) {
        return m_components[i];
    }

    void X(double const c) {
        m_components[0] = c;
    }

    void Y(double const c) {
        m_components[1] = c;
    }

    void Z(double const c) {
        m_components[2] = c;
    }
};


Vec3 operator+(Vec3 const& a, Vec3 const& b) {
    return { a[0] + b[0], a[1] + b[1], a[2] + b[2] };
}

Vec3 operator-(Vec3 const& a, Vec3 const& b) {
    return { a[0] - b[0], a[1] - b[1], a[2] - b[2] };
}

Vec3 operator*(Vec3 const& v, double const c) {
    return { c * v[0], c * v[1], c * v[2] };
}

Vec3 operator*(double const c, Vec3 const& v) {
    return { c * v[0], c * v[1], c * v[2] };
}

Vec3 operator/(Vec3 const& v, double const c) {
    if (c == 0.) {
        return {NAN, NAN, NAN};
    }
    return { v[0] / c, v[1] / c, v[2] / c };
}

double Dot(Vec3 const& a, Vec3 const& b) {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

Vec3 Cross(Vec3 const& a, Vec3 const& b) {
    return { a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0] };
}
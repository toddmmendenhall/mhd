#pragma once

#include <array>

namespace GEOM_UTILS {

class Vector2D {
public:
    Vector2D(double const x, double const y);

private:
    std::array<double, 2> m_components;
};

class Vector3D {
public:
    Vector3D(double const x, double const y, double const z);

    double const X() const;

    double const Y() const;

    double const Z() const;

    void X(double const c);

    void Y(double const c);

    void Z(double const c);

    double const MagnitudeSquared() const;
    
    double const Magnitude() const;

private:
    std::array<double, 3> m_components;
};
    
} // namespace GEOM_UTILS




// Vec3 operator+(Vec3 const& a, Vec3 const& b) {
//     return { a[0] + b[0], a[1] + b[1], a[2] + b[2] };
// }

// Vec3 operator-(Vec3 const& a, Vec3 const& b) {
//     return { a[0] - b[0], a[1] - b[1], a[2] - b[2] };
// }

// Vec3 operator*(Vec3 const& v, double const c) {
//     return { c * v[0], c * v[1], c * v[2] };
// }

// Vec3 operator*(double const c, Vec3 const& v) {
//     return { c * v[0], c * v[1], c * v[2] };
// }

// Vec3 operator/(Vec3 const& v, double const c) {
//     if (c == 0.) {
//         return {NAN, NAN, NAN};
//     }
//     return { v[0] / c, v[1] / c, v[2] / c };
// }

// double Dot(Vec3 const& a, Vec3 const& b) {
//     return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
// }

// Vec3 Cross(Vec3 const& a, Vec3 const& b) {
//     return { a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0] };
// }
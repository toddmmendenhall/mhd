#pragma once

#include <array>

namespace GEOM_UTILS {

class Point2D {
public:
    Point2D(double x, double y);

private:
    std::array<double, 2> m_components;
};

} // namespace GEOM_UTILS

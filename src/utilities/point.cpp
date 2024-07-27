#include "point.hpp"

namespace GEOM_UTILS {

Point2D::Point2D(double x, double y) {
    m_components[0] = x;
    m_components[1] = y;
}

} // namespace GEOM_UTILS

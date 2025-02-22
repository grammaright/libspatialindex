#include <cstring>
#include <cmath>
#include <limits>

#include <spatialindex/SpatialIndex.h>

using namespace SpatialIndex;

double haversineDistance(double lat1, double lon1, double lat2, double lon2) {
  lat1 *= M_PI / 180.0;
  lon1 *= M_PI / 180.0;
  lat2 *= M_PI / 180.0;
  lon2 *= M_PI / 180.0;

  double dlat = lat2 - lat1;
  double dlon = lon2 - lon1;

  double a = std::pow(std::sin(dlat / 2), 2) +
             std::cos(lat1) * std::cos(lat2) * std::pow(std::sin(dlon / 2), 2);
  double c = 2 * std::atan2(std::sqrt(a), std::sqrt(1 - a));

  double radius = 6371.0;  // Earth's radius in kilometers
  return radius * c;
}

double SpherePoint::getMinimumDistance(const IShape& in) const {
  const Point* p = dynamic_cast<const Point*>(&in);
  if (p != nullptr) {
    return getMinimumDistance(*p);
  }

  const Region* r = dynamic_cast<const Region*>(&in);
  if (r == nullptr) {
    throw Tools::IllegalArgumentException(
        "SpherePoint::getMinimumDistance: Unknown shape type.");
  }

  if (m_dimension != 2) {
	throw Tools::IllegalArgumentException(
		"SpherePoint::getMinimumDistance: Only 2D points are supported.");
  }

  // check it is in the box first
  bool is_in_box = false;
  for (int index = 0; index < 2; ++index) {
    auto diff1 = m_pCoords[index] - r->m_pHigh[index];
    auto diff2 = m_pCoords[index] - r->m_pLow[index];

    if (diff1 * diff2 < 0) {
      is_in_box = true;
      return 0;
    }
  }

  // jsut put unrolled loop here
  double dist = std::numeric_limits<double>::infinity();
  dist = std::min(dist, haversineDistance(m_pCoords[1], m_pCoords[0],
                                          r->m_pLow[1], r->m_pLow[0]));
  dist = std::min(dist, haversineDistance(m_pCoords[1], m_pCoords[0],
                                          r->m_pHigh[1], r->m_pLow[0]));
  dist = std::min(dist, haversineDistance(m_pCoords[1], m_pCoords[0],
                                          r->m_pLow[1], r->m_pHigh[0]));
  dist = std::min(dist, haversineDistance(m_pCoords[1], m_pCoords[0],
                                          r->m_pHigh[1], r->m_pHigh[0]));

  return dist;
}

double SpherePoint::getMinimumDistance(const Point& p) const {
  if (m_dimension != p.m_dimension) {
    throw Tools::IllegalArgumentException(
        "Point::getMinimumDistance: Shapes have different number of "
        "dimensions.");
  }

  double lat1 = m_pCoords[1];
  double lon1 = m_pCoords[0];
  double lat2 = p.m_pCoords[1];
  double lon2 = p.m_pCoords[0];

  return haversineDistance(lat1, lon1, lat2, lon2);
}

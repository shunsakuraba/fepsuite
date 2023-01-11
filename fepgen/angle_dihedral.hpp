#pragma once

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <cmath>

inline double dihedral(const Eigen::Vector3d &a1,
		       const Eigen::Vector3d &a2,
		       const Eigen::Vector3d &a3,
		       const Eigen::Vector3d &a4)
{
  Eigen::Vector3d r21 = a1 - a2;
  Eigen::Vector3d r23 = a3 - a2;
  Eigen::Vector3d r34 = a4 - a3;

  r23.normalize();

  Eigen::Vector3d r21v = r21 - r21.dot(r23) * r23;
  Eigen::Vector3d r34v = r34 - r34.dot(r23) * r23;

  r21v.normalize();
  r34v.normalize();
  
  Eigen::Vector3d z = r23.cross(r21v);
  
  double cosphi = r34v.dot(r21v);
  double sinphi = r34v.dot(z);
  
  return atan2(sinphi, cosphi);
}

inline double angle(const Eigen::Vector3d &a1,
		    const Eigen::Vector3d &a2,
		    const Eigen::Vector3d &a3)
{
  Eigen::Vector3d r21 = a1 - a2;
  Eigen::Vector3d r23 = a3 - a2;

  r21.normalize();
  r23.normalize();

  Eigen::Vector3d c = r21.cross(r23);
  
  return atan2(c.norm(), r21.dot(r23));
}


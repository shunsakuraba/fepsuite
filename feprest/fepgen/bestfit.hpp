/*
  Best-fit module, works with Eigen3 or later
 */
#pragma once

/* 
 Copyright (c) 2013, Shun Sakuraba
 All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

3. Neither the name of the Japan Atomic Energy Agency nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <Eigen/Eigenvalues>

static void calc_trans(const Eigen::Matrix3Xd &refstr,
		       const Eigen::Matrix3Xd &targetstr,
		       const Eigen::VectorXd &mass,
		       Eigen::Quaterniond& retq,
		       Eigen::Translation<double, 3>& trans)
{
  int refatoms = refstr.cols();
  int targetatoms = targetstr.cols();
  assert(refatoms == targetatoms || !"Reference and target structure must have the same number of atoms.");
  assert(refstr.rows() == 3 || !"Number of leading dimensions of refstr must be 3");
  assert(targetstr.rows() == 3 || !"Number of leading dimensions of targetstr must be 3");
  assert(mass.rows() == refatoms || !"Number of elements in mass does not match with the number of atoms");

  Eigen::Vector3d targetmean(Eigen::Vector3d::Zero());
  for(int i = 0; i < targetstr.cols(); ++i) {
    targetmean += targetstr.col(i) * mass(i);
  }
  targetmean /= mass.sum();

  Eigen::Matrix3d rtensor = 
    (targetstr.colwise() - targetmean) * mass.asDiagonal() * refstr.transpose();
  
  Eigen::Matrix4d matmax;
  // Eigen documents that it only reads lower half of the matrix, but it isn't (as of 3.1.2)
  matmax(0, 0) = + rtensor(0, 0) + rtensor(1, 1) + rtensor(2, 2);
  matmax(1, 1) = + rtensor(0, 0) - rtensor(1, 1) - rtensor(2, 2);
  matmax(2, 2) = - rtensor(0, 0) + rtensor(1, 1) - rtensor(2, 2);
  matmax(3, 3) = - rtensor(0, 0) - rtensor(1, 1) + rtensor(2, 2);

  matmax(1, 0) = + rtensor(1, 2) - rtensor(2, 1);
  matmax(2, 0) = + rtensor(2, 0) - rtensor(0, 2);
  matmax(3, 0) = + rtensor(0, 1) - rtensor(1, 0);

  matmax(2, 1) = + rtensor(0, 1) + rtensor(1, 0);
  matmax(3, 2) = + rtensor(1, 2) + rtensor(2, 1);
  matmax(3, 1) = + rtensor(2, 0) + rtensor(0, 2);

  matmax(0, 1) = matmax(1, 0);
  matmax(0, 2) = matmax(2, 0);
  matmax(0, 3) = matmax(3, 0);

  matmax(1, 2) = matmax(2, 1);
  matmax(2, 3) = matmax(3, 2);
  matmax(1, 3) = matmax(3, 1);

  Eigen::SelfAdjointEigenSolver<Eigen::Matrix4d> solver(matmax);
  //solver.compute(matmax);
 
  Eigen::Vector4d q = solver.eigenvectors().col(3);
  retq = Eigen::Quaterniond(q(0), q(1), q(2), q(3));
  trans = Eigen::Translation<double, 3>(- targetmean(0), - targetmean(1), - targetmean(2));
}

inline void bestfit(const Eigen::Matrix3Xd &refstr,
	     Eigen::VectorXd mass,
	     Eigen::Matrix3Xd &targetstr)
{
  Eigen::Vector3d refmean(Eigen::Vector3d::Zero());
  for(int i = 0; i < refstr.cols(); ++i) {
    refmean += refstr.col(i) * mass(i);
  }
  refmean /= mass.sum();

  Eigen::Quaterniond q;
  Eigen::Translation<double, 3> t;

  calc_trans((refstr.colwise() - refmean),
	     targetstr,
	     mass,
	     q, t);

  for(int i = 0; i < targetstr.cols(); ++i){
    Eigen::Vector3d c = Eigen::Translation<double, 3>(refmean(0), refmean(1), refmean(2)) * q * t * targetstr.col(i);
    targetstr.col(i) = c;
  }
}

inline void bestfit_all(const Eigen::Matrix3Xd &refstr,
		 const Eigen::VectorXd &mass,
		 Eigen::MatrixXd &targetstr)
{
  Eigen::Vector3d refmean(Eigen::Vector3d::Zero());
  for(int i = 0; i < refstr.cols(); ++i) {
    refmean += refstr.col(i) * mass(i);
  }
  refmean /= mass.sum();

  int ncoord = targetstr.rows();
  assert(ncoord == 3 * refstr.cols() || !"bestfit_all: Number of coordinates in refstr should match");
  int natoms = refstr.cols();

  Eigen::Quaterniond q;
  Eigen::Translation<double, 3> t;

  Eigen::Matrix3Xd refstr0 = (refstr.colwise() - refmean).matrix();

  for(int T = 0; T < targetstr.cols(); ++T) {
    Eigen::Matrix3Xd tsmat = Eigen::Map<Eigen::Matrix3Xd>(targetstr.col(T).data(), 3, natoms);
    calc_trans(refstr0,
	       tsmat,
	       mass,
	       q, t);
    
    for(int i = 0; i < natoms; ++i){
      Eigen::Vector3d c = Eigen::Translation<double, 3>(refmean(0), refmean(1), refmean(2)) * q * t * targetstr.block<3,1>(3 * i, T);
      targetstr.block<3,1>(3 * i, T) = c;
    }
  }
}




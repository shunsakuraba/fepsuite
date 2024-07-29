#include "select.hpp"
#include "bestfit.hpp"

#include <iostream>

using namespace std;
using namespace Eigen;

VectorXd set_selected_mass(const vector<string>& atomnames,
                           const string& sel)
{
  VectorXd ret(VectorXd::Zero(atomnames.size()));
  
  for(int i = 0; i < ret.rows(); ++i) {
    if(atomnames[i] == sel) {
      ret(i) = 1.0;
    }
  }
  
  return ret;
}

void fit_selected(const VectorXd& Amass,
                  const VectorXd& Bmass,
                  const Matrix3Xd& Acoords,
                  Matrix3Xd& Bcoords)
{
  Vector3d refmean(Vector3d::Zero());
  int selcount = 0;
  for(int i = 0; i < Acoords.cols(); ++i) {
    if(Amass(i) > 0.0) {
      refmean += Acoords.col(i) * Amass(i);
      selcount++;
    }
  }
  if(selcount == 0) {
    cerr << "Warning: too few atoms are selected for best-fitting. The program assumes you have best-fitted two PDB files already by some other programs." << endl;
    return;
  }
  refmean /= Amass.sum();
  int selcountB = 0;
  for(int i = 0; i < Bmass.rows(); ++i) {
    if(Bmass(i) > 0) selcountB++;
  }
  if(selcount != selcountB) {
    cerr << "Warning: two PDBs have unequal number of C-alpha or phosphate atoms. The program assumes you have best-fitted two PDB files already by some other programs." << endl;
    // returns without best-fitting
    return;
  }
  if(selcountB == 0) {
    cerr << "Warning: too few atoms are selected for best-fitting. The program assumes you have best-fitted two PDB files already by some other programs." << endl;
    return;
  }
  
  Quaterniond q;
  Translation<double, 3> t;
  
  VectorXd mass(selcount);
  Matrix3Xd refstr(3, selcount);
  Matrix3Xd targetstr(3, selcount);
  int c;
  c = 0;
  for(int i = 0; i < Acoords.cols(); ++i) {
    if(Amass(i) > 0.0) {
      mass(c) = Amass(i);
      refstr.col(c) = Acoords.col(i);
      ++c;
    }
  }
  c = 0;
  for(int i = 0, c = 0; i < Bcoords.cols(); ++i) {
    if(Bmass(i) > 0.0) {
      targetstr.col(c) = Bcoords.col(i);
      ++c;
    }
  }
  
  calc_trans((refstr.colwise() - refmean),
             targetstr,
             mass,
             q, t);
  
  for(int i = 0; i < Bcoords.cols(); ++i){
    Eigen::Vector3d c = Eigen::Translation<double, 3>(refmean(0), refmean(1), refmean(2)) * q * t * Bcoords.col(i);
    Bcoords.col(i) = c;
  }
}

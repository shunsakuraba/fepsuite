#include "assign.hpp"

using namespace std;
using namespace Eigen;

static char atomtype(const string& an)
{
  for(size_t i = 0; i < an.length(); ++i) {
    if(isalpha(an[i])){
      return an[i];
    }
  }
  return 0;
}

void assign_atoms(const string& process_atoms,
                  const vector<string>& Anames,
                  const vector<string>& Bnames,
                  const MatrixXd& distmat,
                  vector<int>& assignBofA,
                  vector<int>& assignAofB,
                  double threshold)
{
  vector<bool> Aenable(Anames.size());
  vector<bool> Benable(Bnames.size());

  for(int i = 0; i < (int)Anames.size(); i++) {
    char at = atomtype(Anames[i]);
    Aenable[i] = at && (process_atoms.find(at) != string::npos);
  }

  for(int i = 0; i < (int)Bnames.size(); i++) {
    char at = atomtype(Bnames[i]);
    Benable[i] = at && (process_atoms.find(at) != string::npos);
  }

  for(int i = 0; i < (int)Anames.size(); ++i) {
    if(!Aenable[i]) continue;
    if(assignBofA[i] != -1) continue;

    int minatm = -1;
    double mindist = numeric_limits<double>::max();
    for(int j = 0; j < (int)Bnames.size(); ++j) {
      if(!Benable[j]) continue;
      if(assignAofB[j] != -1) continue;
      double d = distmat(i, j);
      if(d < mindist){
        minatm = j;
        mindist = d;
      }
    }
    if(minatm >= 0 && mindist < threshold) {
      // found the atom
      assignBofA[i] = minatm;
      assignAofB[minatm] = i;
    }
  }
}

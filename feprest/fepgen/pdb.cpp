#include "pdb.hpp"
#include <fstream>
#include <stdexcept>
#include "string_util.hpp"

using namespace std;
using namespace Eigen;

static
int atoi_str(const string& s)
{
  if(s.length() == 0)
    throw runtime_error("atoi_str: length = 0");

  const char* s_z = s.c_str();
  char *q;
  int retr = strtol(s_z, &q, 10);
  if(*q != '\0')
    throw runtime_error("atoi_str: parse failed");
  return retr;
}

pdb::pdb()
  : altloc('\0')
{
}

pdb::pdb(const string& fname, char altloc)
  : altloc(altloc)
{
  fh.reset(new ifstream(fname.c_str()));

  if(!load_pdb()) 
    throw runtime_error("Failed to read PDB file");
}

bool
pdb::next()
{
  if(!fh) return false;
  return load_pdb();
}

bool
pdb::load_pdb()
{
  ifstream& ifs = *(fh.get());
  
  string line;

  vector<Vector3d> tmpcrds;
  vector<double> tmpbetas;
  vector<double> tmpoccupancies;

  vector<char> tmpchains;
  vector<string> tmpresiduenames;
  vector<int> tmpresids;
  vector<string> tmpatomnames;

  while(getline(ifs, line)) {
    if(line.empty()) continue;
    if(line.size() < 6) continue;
    string cmd = line.substr(0, 6);
    if(cmd == "ATOM  " || cmd == "HETATM") {
      if(line.length() < 54) return false; // Line too short, corrupt

      char altloc_cur = line[16];
      // Skip lines with unmatched altloc value
      if(altloc != '\0' && 
         altloc_cur != ' ' &&
         altloc_cur != altloc)
        continue;

      char chain = line[21];
      string resn = trim(line.substr(17, 3));
      int resid = atoi_str(trim(line.substr(22, 4)));
      string atomn = trim(line.substr(12, 4));
      if(isdigit(atomn[0])) {
        atomn = atomn.substr(1, 3) + atomn.substr(0, 1);
      }

      int r;
      string s = line.substr(30, 24);
      double x, y, z;
      r = sscanf(s.c_str(), "%8lf%8lf%8lf", &x, &y, &z);
      if(r != 3) return false;

      double beta = 0.;
      if(s.length() >= 68) {
        s = line.substr(60, 6);
        r = sscanf(s.c_str(), "%6lf", &beta);
        if(r != 1) return false;
      }

      double occupancy = 0.;
      if(s.length() >= 60) {
        s = line.substr(54, 6);
        r = sscanf(s.c_str(), "%6lf", &occupancy);
        if(r != 1) return false;
      }

      tmpcrds.push_back(Vector3d(x, y, z));
      tmpbetas.push_back(beta);
      tmpoccupancies.push_back(occupancy);

      tmpchains.push_back(chain);
      tmpresiduenames.push_back(resn);
      tmpresids.push_back(resid);
      tmpatomnames.push_back(atomn);
    }else if(cmd == "ENDMDL") {
      break;
    }
  }
  
  size_t cols = tmpcrds.size();

  if(cols == 0) return false;
  
  coords = Matrix3Xd(3, cols);
  betas = VectorXd(cols);
  occupancies = VectorXd(cols);
  
  for(size_t i = 0; i < cols; ++i) {
    coords.col(i) = tmpcrds[i];
    betas(i) = tmpbetas[i];
    occupancies(i) = tmpoccupancies[i];
  }
  
  chains.swap(tmpchains);
  residuenames.swap(tmpresiduenames);
  resids.swap(tmpresids);
  atomnames.swap(tmpatomnames);

  if(!ifs)
    fh.reset();

  return true;
}


#include "pdb.hpp"
#include "cmdline.h"
#include <string>
#include <iostream>
#include <Eigen/Core>
#include <stdexcept>
#include <cassert>
#include <iomanip>
#include <fstream>
#include <stack>

using namespace std;
using namespace Eigen;

const vector<pair<string, double>> vdwtable = {
  {"Cl", 1.75},
  {"Br", 1.85},
  {"H", 1.20},
  {"C", 1.70},
  {"O", 1.52},
  {"N", 1.55},
  {"F", 1.47},
  {"P", 1.80},
  {"S", 1.80},
};


vector<vector<bool>> get_connectivity(const pdb &p)
{
  const Matrix3Xd &coords = p.get_coords();
  VectorXd radii(coords.cols());
  for(int i = 0; i < radii.rows(); ++i) {
    double radius = 1.5; // default
    const string& name = p.get_atomnames()[(size_t)i];
    for(const auto &v: vdwtable) {
      size_t keylen = v.first.length();
      if(name.length() >= keylen &&
         name.substr(0, keylen) == v.first) {
        radius = v.second;
        break;
      }
    }
    radii[i] = radius;
  }

  vector<vector<bool>> ret(coords.cols(), vector<bool>(coords.cols(), false));
  for(int i = 0; i < coords.cols(); ++i) {
    for(int j = 0; j < coords.cols(); ++j) {
      if((coords.col(i) - coords.col(j)).norm() < max(radii[(size_t)i], radii[(size_t)j])) {
        //cerr << "Bond b/w " << i << " " << j << endl; 
        ret[i][j] = true;
      }
    }
  }
  return ret;
}

bool match_impl(size_t curptr, vector<size_t> &input2ref,
                vector<bool> &visited,
                const vector<size_t> &inp_reorder,
                const vector<vector<bool>>& refconn,
                const vector<vector<bool>>& inpconn,
                const vector<vector<bool>>& match_atomtype)
{
  if(curptr >= inpconn.size()) {
    return true;
  }

  assert(input2ref.size() == curptr);

  for(size_t trial = 0; trial < refconn.size(); trial++){
    // try to set input2ref[curptr] to be trial
    if(visited[trial]) continue;
    bool isok = true;
    for(size_t prev = 0; prev < curptr; prev++) {
      size_t atom = input2ref[prev];
      if(inpconn[inp_reorder[curptr]][inp_reorder[prev]] != refconn[trial][atom]){
        isok = false;
        break;
      }
    }
    if(!isok) continue;
    
    // actually set and try search
    input2ref.push_back(trial);
    visited[trial] = true;
    isok = match_impl(curptr + 1, input2ref, visited, inp_reorder, refconn, inpconn, match_atomtype);
    if(isok) {
      return true;
    }
    visited[trial] = false;
    input2ref.pop_back();
  }
  return false;
}

vector<size_t> match_atoms(const pdb &ref, const pdb &inp)
{
  vector<size_t> input2ref;

  size_t nref = ref.get_atomnames().size();
  size_t ninp = inp.get_atomnames().size();
  vector<vector<bool>> refconn = get_connectivity(ref);
  vector<vector<bool>> targetconn = get_connectivity(inp);

  vector<vector<bool>> match_atomtype(nref, vector<bool>(ninp, false));
  for(size_t i = 0; i < nref; ++i) {
    for(size_t j = 0; j < ninp; ++j) {
      bool match;
      const string& refname = ref.get_atomnames()[i];
      const string& inpname = inp.get_atomnames()[j];
      if(refname.length() >= 2 && refname.substr(0, 2) == "Cl") {
        match = inpname.length() >= 2 && inpname.substr(0, 2) == "Cl";
      }else{
        match = (refname[0] == inpname[0]);
      }
      match_atomtype[i][j] = match;
    }
  }

  vector<size_t> degree(ninp, 0);
  for(size_t i = 0; i < ninp; ++i) {
    int count = 0;
    for(size_t j = 0; j < ninp; ++j) {
      count += (int) targetconn[i][j];
    }
    degree[i] = count - 1; // -1 for self edge
  }

  size_t smallest_degree = 999;
  size_t initialnode = 0;
  for(size_t i = 0; i < ninp; ++i) {
    if(degree[i] < smallest_degree) {
      smallest_degree = degree[i];
      initialnode = i;
    }
  }

  // DFS to make node reordering
  vector<size_t> inp_reorder;
  {
    vector<bool> visited(ninp, false);
    stack<size_t> dfsstack;
    dfsstack.push(initialnode);

    while(!dfsstack.empty()) {
      size_t t = dfsstack.top();
      dfsstack.pop();
      if(visited[t]) continue;
      visited[t] = true;
      inp_reorder.push_back(t);
      for(size_t j = 0; j < ninp; j++) {
        if(visited[j]) continue;
        if(!targetconn[t][j]) continue;
        dfsstack.push(j);
      }
    }
    assert(inp_reorder.size() == ninp);
  }

  // TODO pre-sort by connectivity to get faster assignment
  vector<bool> visited(ref.get_atomnames().size(), false);
  bool found = match_impl(0, input2ref, visited, inp_reorder, refconn, targetconn, match_atomtype);
  
  if(!found) {
    throw runtime_error("No solution found");
  }

  vector<size_t> input2ref_real(ninp);
  for(size_t i = 0; i < ninp; ++i) {
    input2ref_real[inp_reorder[i]] = input2ref[i];
  }
  return input2ref_real;
}

int main(int argc, char* argv[])
{
  cmdline::parser p;

  p.add("help", 0, "Print this message");
  p.add<string>("reference", 'r', "Reference PDB", true);
  p.add<string>("input", 'f', "Input PDB", true);
  p.add<string>("output", 'o', "Output PDB", true);
  p.add("reorder-to-reference", 0, "Reorder to refernce structure (input and reference must have same number of atoms");

  {
    bool ok = p.parse(argc, argv);
    if(!ok || p.exist("help")){
      cerr << "Match two atom names in two PDB files" << endl;
      cerr << p.usage();
      return ((ok && p.exist("help")) ? EXIT_SUCCESS : EXIT_FAILURE);
    }
  }

  pdb refpdb(p.get<string>("reference"));
  pdb targetpdb(p.get<string>("input"));

  ofstream crdfs(p.get<string>("output"));
  crdfs.setf(ios::fixed);
  vector<size_t> input2ref = match_atoms(refpdb, targetpdb);
  vector<size_t> printorder(input2ref.size());
  if(p.exist("reorder-to-reference")) {
    if(refpdb.get_atomnames().size() != targetpdb.get_atomnames().size()) {
      cerr << "To reorder two PDBs must have matched number of atoms" << endl;
      return EXIT_FAILURE;
    }
    for(size_t i = 0; i < input2ref.size(); ++i) {
      printorder[input2ref[i]] = i;
    }
  }else{
    for(size_t i = 0; i < input2ref.size(); ++i) {
      printorder[i] = i;
    }
  }
  size_t atomcount = 1;
  for(size_t i: printorder) {
    size_t ref = input2ref[i];
    crdfs << setw(6) << "ATOM  " 
          << setw(5) << std::right << atomcount++
          << " ";
    const string& atom = refpdb.get_atomnames()[ref];
    if(atom.length() == 4) {
      if(isdigit(atom[3])) {
        crdfs << atom.substr(3, 1) << setw(3) << std::left
              << atom.substr(0, 3);
      }else{
        crdfs << atom;
      }
    }else{
      crdfs << " "
            << setw(3) << std::left
            << atom.substr(0,3);
    }
    const string& resname = refpdb.get_residuenames()[ref];
    int resid = refpdb.get_resids()[ref];
    Vector3d crd = targetpdb.get_coords().col(i);
    crdfs << " " // altLoc
          << setw(3) << std::right << resname
          << " "
          << "A" // chainID
          << setw(4) << std::right << resid
          << " " // icode
          << "   "
          << setw(8) << std::right << setprecision(3) << crd(0)
          << setw(8) << std::right << setprecision(3) << crd(1)
          << setw(8) << std::right << setprecision(3) << crd(2)
          << endl;
  }
  crdfs << "TER   " << endl;
}

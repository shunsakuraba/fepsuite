#include <string>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <cmath>
#include <algorithm>
#include "topology.hpp"

#include <iostream> // debug

using namespace std;

topology::topology()
{
}

topology::topology(const string& fname)
  : defaults(AMBER), nexcl(0)
{
  ifstream ifs(fname.c_str());
  
  string line;
  string state;
  int nmoltype = 0;
  while(getline(ifs, line)) {
    if(line == ""){
      continue;
    }

    if(line[0] == '[') {
      istringstream is(line);
      char bopen, bclose;
      string cmd;
      is >> bopen >> cmd >> bclose;
      if(is.bad() || bopen != '[' || bclose != ']') {
        throw runtime_error("topology::topology bracket parse error");
      }

      state = cmd;
      continue;
    }
    if(line[0] == '#') {
      // currently ignore includes / defines ...
      continue;
    }

    {
      size_t p = line.find(';');
      if(p != string::npos) {
        line = line.substr(0, p);
      }
    }
    {
      size_t p = line.find_last_not_of(" \t");
      if(p != string::npos) {
        line = line.substr(0, p + 1);
      }
    }
    if(line == ""){
      continue;
    }

    if(state == "defaults") {
      istringstream is(line);
      int nbfunc, comb_rule;
      string gen_pairs;
      double fudgeLJ, fudgeQQ;
      is >> nbfunc >> comb_rule >> gen_pairs >> fudgeLJ >> fudgeQQ;
      if(nbfunc == 1 && 
         comb_rule == 2 &&
         gen_pairs == "yes" &&
         fabs(fudgeLJ - 0.5) < 1e-3 &&
         fabs(fudgeQQ - 0.8333) < 1e-3) {
        defaults = AMBER;
      }else if(nbfunc == 1 &&
               comb_rule == 2 &&
               gen_pairs == "yes" &&
               fabs(fudgeLJ - 1.0) < 1e-3 &&
               fabs(fudgeQQ - 1.0) < 1e-3) {
        defaults = CHARMM;
      }else if(nbfunc == 1 &&
               comb_rule == 3 &&
               gen_pairs == "yes" &&
               fabs(fudgeLJ - 0.5) < 1e-3 &&
               fabs(fudgeQQ - 0.5) < 1e-3) {
        defaults = OPLS;
      }else{
        throw runtime_error("topology::topology: unknown default type");
      }
    }else if(state == "atomtypes") {
      istringstream is(line);
      string atype;
      is >> atype;
      atomtypes[atype] = line;
    }else if(state == "pairtypes") {
      istringstream is(line);
      string atype, btype;
      int func;
      is >> atype >> btype >> func;
      pairtypes[make_tuple(atype, btype, func)] = line;
    }else if(state == "moleculetype") {
      istringstream is(line);
      is >> moleculename >> nexcl;
      if(is.fail() || is.bad()) {
        throw runtime_error("moleculetype format error");
      }
      ++nmoltype;
      if(nmoltype != 1) {
        throw runtime_error("Unsupported: moleculetype appeared more than once");
      }
    }else if(state == "atoms") {
      istringstream is(line);
      int nr, resi, cgnr;
      string type, res, atom;
      double charge, mass;
      is >> nr >> type >> resi >> res >> atom >> cgnr >> charge >> mass;
      if(is.fail() || is.bad()) {
        throw runtime_error("atoms section format error");
      }
      names.push_back(atom);
      types.push_back(type);
      resnames.push_back(res);
      resids.push_back(resi);
      charges.push_back(charge);
      masses.push_back(mass);
    }else if(state == "pairs") {
      istringstream is(line);
      int a, b, func;
      is >> a >> b >> func;
      if(is.fail() || is.bad()) {
        throw runtime_error("pairs format error");
      }
      pairs[make_pair(a - 1, b - 1)] = func;
    }else if(state == "bonds") {
      istringstream is(line);
      int a, b, func;
      is >> a >> b >> func;
      if(is.fail() || is.bad()) {
        throw runtime_error("bonds format error");
      }
      vector<double> &bond = bonds[make_tuple(a - 1, b - 1, func)];
      
      double v;
      while(is >> v){
        bond.push_back(v);
      }
    }else if(state == "angles") {
      istringstream is(line);
      int a, b, c, func;
      is >> a >> b >> c >> func;
      if(is.fail() || is.bad()) {
        throw runtime_error("angles format error");
      }
      vector<double> &angle = angles[make_tuple(a - 1, b - 1, c - 1, func)];
      
      double v;
      while(is >> v){
        angle.push_back(v);
      }
    }else if(state == "dihedrals") {
      istringstream is(line);
      int a, b, c, d, func;
      is >> a >> b >> c >> d >> func;
      if(is.fail() || is.bad()) {
        throw runtime_error("dihedrals format error");
      }
      vector<double> vals;
      double v;
      while(is >> v){
        vals.push_back(v);
      }
      int addendum = 0;
      if(func == 1) {
          addendum = (int)vals.back(); // rep values
      }
      if(addendum == 0 && diheds.count(make_tuple(a - 1, b - 1, c - 1, d - 1, func, addendum)) > 0) {
        throw runtime_error("unsupported: dihedral multiple entries");
      }
      vector<double> &dihedral = diheds[make_tuple(a - 1, b - 1, c - 1, d - 1, func, addendum)];
      dihedral = vals;
    }else if(state == "system" || state == "molecules") {
      // do nothing
    }else{
      throw runtime_error((string("Unsupported section ") + state).c_str());
    }
  }

  if(ifs.bad()) {
    throw runtime_error("I/O error on topology constructor");
  }
}

void
topology::write(const string& fname)
{
  ofstream ofs(fname.c_str());

  if(ofs.bad()) {
    throw runtime_error("I/O error on topology::write");
  }
  throw runtime_error("Unimplemented: toplogy::write()");
}

void
topology::convert_bonds_to_adj_list(vector<vector<int> > &adj_list) const
{
  adj_list.clear();
  adj_list.resize(this->names.size());
  for(const auto& bonditer: this->bonds) {
    int a = get<0>(bonditer.first);
    int b = get<1>(bonditer.first);
    adj_list[a].push_back(b);
    adj_list[b].push_back(a);
  }
  for(auto &v: adj_list) {
    sort(v.begin(), v.end());
  }
}


void
topology::convert_pairs_to_adj_list(vector<vector<int> > &list) const
{
  list.clear();
  list.resize(this->names.size());
  for(const auto& pairiter: this->pairs) {
    int a = pairiter.first.first;
    int b = pairiter.first.second;
    list[a].push_back(b);
    list[b].push_back(a);
  }
  for(auto &v: list) {
    sort(v.begin(), v.end());
  }
}



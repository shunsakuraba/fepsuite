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
      cerr << "Warning: topology file is not preprocessed!" << endl;
      continue;
    }
    if(line[0] == '*') {
      // another form of comment.
      continue;
    }

    {
      size_t p = line.find(';');
      if(p != string::npos) {
        line = line.substr(0, p);
      }
    }
    {
      size_t p = line.find_last_not_of(" \t\f\v\n\r");
      if(p != string::npos) {
        line.erase(p + 1);
      }else{
        line.clear(); // all whitespace
      }
    }
    if(line == ""){
      continue;
    }

    bool gen_pairs;
    if(state == "defaults") {
      istringstream is(line);
      int nbfunc, comb_rule;
      string gen_pairs_str;
      double fudgeLJ, fudgeQQ;
      is >> nbfunc >> comb_rule >> gen_pairs_str >> fudgeLJ >> fudgeQQ;
      if(!is){
          // XXX: fudgeLJ and fudgeQQ are optional
          throw runtime_error("Failed to parse [ defaults ] section");
      }
      if(gen_pairs_str == "yes"){
          gen_pairs = true;
      }

      if(nbfunc == 1 && 
         comb_rule == 2 &&
         gen_pairs &&
         fabs(fudgeLJ - 0.5) < 1e-3 &&
         fabs(fudgeQQ - 0.8333) < 1e-3) {
        defaults = AMBER;
      }else if(nbfunc == 1 &&
               comb_rule == 2 &&
               gen_pairs &&
               fabs(fudgeLJ - 1.0) < 1e-3 &&
               fabs(fudgeQQ - 1.0) < 1e-3) {
        defaults = CHARMM;
      }else if(nbfunc == 1 &&
               comb_rule == 3 &&
               gen_pairs &&
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

      // check for badly formatted ones (typically amber14sb.ff)
      if(isdigit(atype[0])) {
        cerr << "[ atomtypes ] contains type " << atype << " which starts from non-alphabets." << endl;
        cerr << "This is invalid atomtype in GROMACS (but some user-generated force field, such as amber14sb.ff, contains the entries because AmberTools allow such kind of atom types.)" << endl;
        cerr << "Consider renaming all such atomtypes to start from alphabets, e.g. 2C -> A2C." << endl;
        throw runtime_error("topology::topology: invalid force field file detected");
        // Do you know why I added this error message? Otherwise I will be messed with such ill-formed force fields.
      }

      // Old-type [ atomtypes ] may still exist, just to be sure
      vector<string> inputs;
      while(is){
        string token;
        is >> token;
        if(is) inputs.push_back(token);
      }
      if(inputs.size() < 6) {
        cerr << "[ atomtypes ] line is too short" << endl;
        cerr << "\"" << line << "\"" << endl;
        throw runtime_error("topology::topology: failed to parse");
      }
      bool have_bonded_type, have_atomic_number;
      if(inputs[5].length() == 1 && isalpha(inputs[5][0])) {
        have_bonded_type = true;
        have_atomic_number = true;
      }else if(inputs[3].length() == 1 && isalpha(inputs[3][0])) {
        have_bonded_type = false;
        have_atomic_number = false;
      }else{
        have_bonded_type = isalpha(inputs[1][0]);
        have_atomic_number = ! have_bonded_type;
        if(!have_bonded_type) {
          // atomic number should be integer
          if(!std::all_of(inputs[1].begin(), inputs[1].end(), [](char c) { return isdigit(c); })) {
            cerr << "In [ atomtypes ], type " << atype << ": " << endl;
            cerr << "The second entry starts from numerics. This is not a valid atom type for bonds. " << endl;
            throw runtime_error("topology::topology: bondtype assignment error");
          }
        }
      }

      (void) have_atomic_number; // set but unused

      string bond_type = atype;
      if(have_bonded_type) bond_type = inputs[1];
      
      bondatomtypes[atype] = bond_type;
    }else if(state == "pairtypes") {
      istringstream is(line);
      string atype, btype;
      int func;
      is >> atype >> btype >> func;
      pairtypes[make_tuple(atype, btype, func)] = line;
    }else if(state == "bondtypes") {
      istringstream is(line);
      string atype, btype;
      int func;
      is >> atype >> btype >> func;
      bondtypes[make_tuple(atype, btype, func)].push_back(line);
    }else if(state == "angletypes") {
      istringstream is(line);
      string atype, btype, ctype;
      int func;
      is >> atype >> btype >> ctype >> func;
      angletypes[make_tuple(atype, btype, ctype, func)].push_back(line);
    }else if(state == "dihedraltypes") {
      istringstream is(line);
      string atype, btype, ctype, dtype;
      int func;
      is >> atype >> btype >> ctype >> dtype >> func;
      dihedraltypes[make_tuple(atype, btype, ctype, dtype, func)].push_back(line);
    }else if(state == "cmaptypes") {
      // this is mostly unnecessary, but to 0-filld cmaptype if one of state do not have cmap potential.
      istringstream is(line);
      string atype, btype, ctype, dtype, etype;
      int func;
      is >> atype >> btype >> ctype >> dtype >> etype >> func;
      // cmaptypes have strange format using backslash at the end of line
      vector<double> values;
      while(true){
        double v;
        is >> v;
        if(is.fail()) {
          // test whether we can read one more character and it is '\'
          is.clear();
          char c;
          is >> c;
          if(is.fail()) {
            // no characters, i.e., end of the block
            break;
          }
          if(c == '\\') {
            getline(ifs, line);
            is.str(line);
            continue;
          }else{
            throw runtime_error((std::string("Unexpected [ cmaptypes ] character: ") + c).c_str());
          }
        }
        values.push_back(v);
      }
      if(values.size() < 2) {
        throw runtime_error("Invalid [ cmaptypes ] entry: size not found");
      }
      int phi_split = int(values[0]);
      int psi_split = int(values[1]);
      values.erase(values.begin());
      values.erase(values.begin());
      if(values.size() != (size_t)(phi_split * psi_split)) {
        throw runtime_error("Invalid [ cmaptypes ] entry: too few arguments");        
      }
      topology::cmaptype key = std::make_tuple(atype, btype, ctype, dtype, etype, func);
      cmaptypes[key] = std::make_tuple(phi_split, psi_split, values);
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
      if(func == 1 || func == 9) {
        if(vals.empty()) {
          throw runtime_error("unsupported: input must be canonicalized");
        }
        addendum = (int)vals.back(); // rep values
      }
      if(addendum == 0 && diheds.count(make_tuple(a - 1, b - 1, c - 1, d - 1, func, addendum)) > 0) {
        throw runtime_error("unsupported: dihedral multiple entries");
      }
      vector<double> &dihedral = diheds[make_tuple(a - 1, b - 1, c - 1, d - 1, func, addendum)];
      dihedral = vals;
    }else if(state == "cmap") {
      istringstream is(line);
      int a, b, c, d, e, func;
      is >> a >> b >> c >> d >> e >> func;

      cmaps.insert(make_tuple(a - 1, b - 1, c - 1, d - 1, e - 1, func));
    }else if(state == "system" || state == "molecules") {
      // do nothing
    }else if(state == "constrainttypes") {
      // do nothing. constrainttypes are not used unless [ constraints ] section appears, so safely ignored.
    }else if(state == "implicit_genborn_params") {
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



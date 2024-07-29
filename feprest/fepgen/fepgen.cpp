#include "cmdline.h"
#include "bestfit.hpp"
#include "pdb.hpp"
#include "topology.hpp"
#include "assign.hpp"
#include "select.hpp"
#include "angle_dihedral.hpp"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <limits>
#include <set>
#include <queue>
#include <algorithm>
#include <array>
#include "disjoint_set.hpp"

using namespace std;
using namespace Eigen;

enum {
  HAM_CHARGE = 0,
  HAM_VDW,
  HAM_BONDANGLE,
  HAM_DIHEDRAL,
  HAM_DUMMY_RESTRAINT,
  HAM_NR
};

struct dummy_anchoring_parameters
{
  double force_bond;
  double force_angular;
  double force_dihedral;
};

bool verbose = false;
bool quiet = false;

const double nm_of_angstrom = 0.1;

/*
  Returns all indices that has distance exactly n.
  Note that for the case of 5-membered ring with n=3 there are none.
  Returned array is indexed by 'A' notation, so you might have to lift it to 'O' by AofO / OofA.
 */
vector<int> list_n_step_adjacency(const vector<vector<int> >& adj,
                                  int n,
                                  int atomix)
{
  vector<int> ret;
  vector<bool> pushed(adj.size(), false);

  queue<pair<int, int> > bfs_queue;

  pushed[atomix] = true;
  bfs_queue.push(make_pair(atomix, 0));

  while(!bfs_queue.empty()) {
    pair<int, int> e = bfs_queue.front();
    bfs_queue.pop();
    int v = e.first;
    int dist = e.second;
    if(dist == n) {
      ret.push_back(v);
      continue;
    }
    for(int x: adj[v]) {
      if(pushed[x]) continue;
      int newdist = dist + 1;
      bfs_queue.push(make_pair(x, newdist));
      pushed[x] = true;
    }
  }

  return ret;
}

/*
  Calculate the integer distance to non-phantom atoms for each atom.
  [top] is the topology and [phantoms] are indices of phantom atoms that do not exist in the /other/ state.
  i.e., if top is Atop, [phantoms] are atoms existing in A only and not in state B.
  Returned array is indexed by 'A' notation, so you might have to lift it to 'O' by AofO / OofA.
 */
vector<int> calculate_int_distance_nonphantom(const topology& top,
                                              const vector<int>& phantoms)
{
  int large = 9999;
  vector<int> distances(top.names.size(), 0);
  for(int p: phantoms){
    distances[p] = large;
  }
  
  vector<vector<int>> adj;
  top.convert_bonds_to_adj_list(adj);

  // um, I know this is a bad implementation but it should not be[TM] so bad...
  // FIXME TODO: rewrite as a BFS from non-phantoms
  queue<int> to_be_processed;
  for(int p: phantoms) {
    to_be_processed.push(p);
  }

  while(!to_be_processed.empty()) {
    int p = to_be_processed.front();
    to_be_processed.pop();
    int lowest_level = large;
    for(int v: adj[p]) {
      if(distances[v] < lowest_level) {
        lowest_level = distances[v];
      }
    }
    if(lowest_level == large) {
      to_be_processed.push(p); // repush to the end
    }else{
      distances[p] = lowest_level + 1;
    }
  }
  
  return distances;
}

static
array<int, 3> find_zmat_like_connectivity_impl(const vector<vector<int> > &adj,
                                               const vector<int>& phantoms,
                                               const vector<int>& int_distances,
                                               const vector<string>& atomnames,
                                               bool allow_hydrogen,
                                               int atomix,
                                               bool *found)
{
  array<int, 3> ret = {0, 0, 0};
  int curdist = int_distances[atomix];
  for(int b: adj[atomix]) {
    if(b == atomix || !(int_distances[b] == 0 || int_distances[b] < curdist)) {
      continue;
    }
    ret[0] = b;
    for(int a: adj[b]) {
      if(a == atomix || a == b || !(int_distances[a] == 0 || int_distances[a] < int_distances[b])) {
        continue;
      }
      ret[1] = a;
      for(int d: adj[a]) {
        if(d == atomix || d == a || d == b ||
           !(int_distances[d] == 0 || int_distances[d] < int_distances[a])) {
          continue;
        }
        if((!allow_hydrogen) && (atomnames[d][0] == 'H')) {
          continue;
        }
        ret[2] = d;
        *found = true;
        return ret;
      }
    }
  }
  // Not found
  *found = false;
  return ret;
}

static
array<int, 3> find_zmat_like_connectivity(const vector<vector<int> > &adj,
                                          const vector<int>& phantoms,
                                          const vector<int>& int_distances,
                                          const vector<string>& atomnames,
                                          int atomix)
{
  bool found;
  array<int, 3> ret;
  ret = find_zmat_like_connectivity_impl(adj, phantoms, int_distances, atomnames, false, atomix, &found);
  if(found) {
    return ret;
  }else{
    // retry allowing H
    ret = find_zmat_like_connectivity_impl(adj, phantoms, int_distances, atomnames, true, atomix, &found);
    if(!found) {
      throw runtime_error("find_zmat_like_connectivity: failed to find");
    }
    return ret;
  }
}

                   

template<typename F>
int find_closest(const Matrix3Xd &Acoords, 
                 const vector<int>& assignAofO,
                 int x,
                 const F& filter)
{
  int N = (int)assignAofO.size();
  int curmin = -1;
  double distmin = numeric_limits<double>::max();
  Vector3d xpt = Acoords.col(assignAofO[x]);
  for(int y = 0; y < N; ++y) {
    if(!filter(y)) {
      continue;
    }
    Vector3d pt = Acoords.col(assignAofO[y]);
    double d = (xpt - pt).norm();
    
    if(d < distmin) {
      distmin = d;
      curmin = y;
    }
  }
  return curmin;
}

static const vector<double> get_bond_entry(const topology &top, topology::bondkeytype& key)
{
  assert(top.bonds.count(key) > 0);
  const vector<double> &entry = top.bonds.at(key);
  if(!entry.empty()) {
    return entry;
  }else{
    // read from bondtype
    const string& at1 = top.types[std::get<0>(key)];
    const string& at2 = top.types[std::get<1>(key)];
    const string& bat1 = top.bondatomtypes.at(at1);
    const string& bat2 = top.bondatomtypes.at(at2);
    topology::bondtype btk(bat1, bat2, std::get<2>(key));
    const vector<string>& ret_str = top.bondtypes.at(btk);
    // FIXME we need to store bondtypes by vector<double>, otherwise we have too many memory assignments here
    vector<double> ret;
    for(const string& v: ret_str) {
      ret.push_back(stod(v));
    }
    return ret;
  }
}

/*
Returns expected bond base distance between two atoms
*/
static double topological_min_dist(const topology &Atop, int a, int b)
{
  for(int functype: {1, 2, 3, 4, 6}) { // bond, G96bond, Morse, cubic, harmonic restraint
    for(int order: {0, 1}) {
      topology::bondkeytype key;
      if(order == 0) {
        key = topology::bondkeytype(a, b, functype);
      }else{
        key = topology::bondkeytype(b, a, functype);
      }
      if(Atop.bonds.count(key) > 0) {
        const vector<double> entry = get_bond_entry(Atop, key);
        assert(entry.size() > 0);
        return entry[0] * 10.0; // for angle types 1-4 & 6 first element is b0 (nm), returns in angstrom
      }
    }
  }
  cerr << "Error: failed to find distance between" << Atop.names[a] << " - " << Atop.names[b] << endl;
  exit(1);
}

static const vector<double> get_angle_entry(const topology &top, topology::anglekeytype& key)
{
  assert(top.angles.count(key) > 0);
  const vector<double> &entry = top.angles.at(key);
  if(!entry.empty()) {
    return entry;
  }else{
    // read from angletype
    const string& at1 = top.types[std::get<0>(key)];
    const string& at2 = top.types[std::get<1>(key)];
    const string& at3 = top.types[std::get<2>(key)];
    const string& bat1 = top.bondatomtypes.at(at1);
    const string& bat2 = top.bondatomtypes.at(at2);
    const string& bat3 = top.bondatomtypes.at(at3);
    topology::angletype btk(bat1, bat2, bat3, std::get<3>(key));
    const vector<string>& ret_str = top.angletypes.at(btk);
    // FIXME we need to store angletypes by vector<double>, otherwise we have too many memory assignments here
    vector<double> ret;
    for(const string& v: ret_str) {
      ret.push_back(stod(v));
    }
    return ret;
  }
}

/*
Returns expected bond base angle between two atoms
*/
static double topological_min_angle(const topology &Atop, int a, int b, int c)
{
  for(int functype: {1, 2, 5, 10 }) { // harmonic, g96, U-B, restraint
    for(int order: {0, 1}) {
      topology::anglekeytype key;
      if(order == 0) {
        key = topology::anglekeytype(a, b, c, functype);
      }else{
        key = topology::anglekeytype(c, b, a, functype);
      }
      if(Atop.angles.count(key) > 0) {
        const vector<double> entry = get_angle_entry(Atop, key);
        assert(entry.size() > 0);
        return entry[0]; // for angle types 1,2,5 and 10 first element is theta0 (degree)
      }
    }
  }
  cerr << "Error: failed to find angle between" << Atop.names[a] << " - " << Atop.names[b] << " - " << Atop.names[c] << endl;
  exit(1);
}

static const vector<double> get_dihedral_entry(const topology &top, topology::dihedkeytype& key)
{
  assert(top.diheds.count(key) > 0);
  const vector<double> &entry = top.diheds.at(key);
  if(!entry.empty()) {
    return entry;
  }else{
    // read from angletype
    const string& at1 = top.types[std::get<0>(key)];
    const string& at2 = top.types[std::get<1>(key)];
    const string& at3 = top.types[std::get<2>(key)];
    const string& at4 = top.types[std::get<3>(key)];
    const string& bat1 = top.bondatomtypes.at(at1);
    const string& bat2 = top.bondatomtypes.at(at2);
    const string& bat3 = top.bondatomtypes.at(at3);
    const string& bat4 = top.bondatomtypes.at(at4);
    topology::dihedraltype dtk(bat1, bat2, bat3, bat4, std::get<3>(key));
    const vector<string>& ret_str = top.dihedraltypes.at(dtk);
    // FIXME we need to store angletypes by vector<double>, otherwise we have too many memory assignments here
    vector<double> ret;
    for(const string& v: ret_str) {
      ret.push_back(stod(v));
    }
    return ret;
  }
}

/*
FIXME untested
 */
static vector<double> scan_rb(const vector<double>& rbfactor)
{
  int slice = 36;
  vector<double> pots;
  for(int i = 0; i < slice; ++i) {
    double angle = 2.0 * M_PI * i / slice;
    double value = 0;
    double cosphi = cos(angle);
    for(int degree = 0; degree < 6; ++degree) {
      value += rbfactor[degree] * pow(cosphi, degree);
    }
    pots.push_back(value);
  } 
  double thresh = *std::min_element(pots.begin(), pots.end()) + 5e-1; // with this 0.5kJ threshold it only accepts close to global minimum
  pots.push_back(pots[0]);
  pots.push_back(pots[1]);
  vector<double> rets;
  for(int i = 1; i <= slice; ++i) {
    if(pots[i] < thresh && pots[i - 1] >= pots[i] && pots[i] <= pots[i + 1]) {
      rets.push_back((i % slice) * 360.0 / slice);
    }
  }
  return rets;
}

/*
FIXME untested
 */
static vector<double> scan_mult(const vector<vector<double> >& entry)
{
  int slice = 36;
  vector<double> pots;
  for(int i = 0; i < slice; ++i) {
    double angle = 2.0 * M_PI * i / slice;
    double value = 0;
    for(const vector<double> &e: entry) {
      double phi = e[0] * M_PI / 180.0; // radian
      double k = e[1];
      int mult = (int) e[2];
      // k (1 + cos(n phi - phi_s))
      value += k * (1. + cos(mult * angle - phi));
    }
    pots.push_back(value);
  }
  double thresh = *std::min_element(pots.begin(), pots.end()) + 5e-1; // with this 0.5kJ threshold it only accepts close to global minimum
  pots.push_back(pots[0]);
  pots.push_back(pots[1]);
  vector<double> rets;
  for(int i = 1; i <= slice; ++i) {
    if(pots[i] < thresh && pots[i - 1] >= pots[i] && pots[i] <= pots[i + 1]) {
      rets.push_back((i % slice) * 360.0 / slice);
    }
  }
  return rets;
}


/*
Returns expected bond base angle between two atoms
*/
static vector<double> topological_min_dihed(const topology &Atop, int a, int b, int c, int d)
{
  static bool warned = false;
  const int mult_max_degree = 6; // scans up to n=6
  for(int functype: {1, 2, 3, 4, 5, 9, 10}) { // supports proper, improper, R-B, periodic improper, fourier, mult. proper, restricted
    for(int order: {0, 1}) {
      vector<vector<double>> found;
      int multmax = (functype == 9 ? mult_max_degree : 0);
      int multmin = (functype == 9 ? 1 : 0);
      for(int addenda = multmin; addenda <= multmax; ++addenda) {
        topology::dihedkeytype key;
        if(order == 0) {
          key = topology::dihedkeytype(a, b, c, d, functype, addenda);
        }else{
          key = topology::dihedkeytype(d, c, b, a, functype, addenda);
        }
        if(Atop.diheds.count(key) > 0) {
          const vector<double> entry = get_dihedral_entry(Atop, key);
          assert(entry.size() > 0);
          found.emplace_back(entry); // for angle types 1,2,5 and 10 first element is theta0
        }
      }
      if(!found.empty()) {
        // return best angle
        vector<double> ret;
        switch(functype) {
          case 1: // proper
          case 4: // periodic improper
            // argmin phi [1 + cos (n phi - phi_s)]
            // <=> n phi - phi_s == 180 (mod 360) # degrees
            // <=> n phi = 180 + phi_s + 360 m
            // <=> phi = (360m + 180 + phi_s) / n
            assert(found.size() == 1);
            {
              double phi = found[0][0];
              int mult = found[0][2]; // phi k mult
              for(int m = 0; m < mult; ++m) {
                ret.push_back((360.0 * m + 180.0 + phi) / mult);
              }
            }
            return ret;
            break;
          case 2: // improper
          case 10: // restraint
            assert(found.size() == 1);
            ret.push_back(found[0][0]);
            return ret;
            break;
          case 3: // R-B
          case 5: // Fourier
            assert(found.size() == 1);
            {
              vector<double> RBfactor;
              if(functype == 3){
                RBfactor = found[0];
              }else{
                const vector<double>& ff = found[0];
                RBfactor = {
                  ff[2] + 0.5 * (ff[1] + ff[3]),
                  0.5 * (-ff[1] + 3 * ff[3]),
                  -ff[2] + 4 * ff[4],
                  -2 * ff[3],
                  -4 * ff[4],
                  0
                }; // gromacs manual 
              }
              return scan_rb(RBfactor);
            }
            break;
          case 9: // multiple
            return scan_mult(found);
            break;
          default:
            if(!warned) {
              cerr << "Warninig: unsupported dihedral type called in topological_min_dihed" << endl;
              warned = true;
            }
        }
      }
    }
  }
  cerr << "Error: failed to find angle between" << Atop.names[a] << " - " << Atop.names[b] << " - " << Atop.names[c] << " " << Atop.names[d] << endl;
  exit(1);
}

// Returns true if angle is defined and close to 180 degree. 
bool unfavored_for_dihedral_restraint(const Matrix3Xd &coords, int a1, int a2, int a3)
{
  if(a1 == -1 || a2 == -1 || a3 == -1) {
    return false;
  }
  const double linear_angle_threshold = 170.0 * M_PI / 180.0;

  double a = angle(coords.col(a1), coords.col(a2), coords.col(a3));
  return (a > linear_angle_threshold);
}

void output_dummies(ofstream &Ofs, const Matrix3Xd &Acoords,
                    double force_bond,
                    double force_angluar,
                    double force_dihedral,
                    const vector<int>& assignAofO,
                    const vector<int>& assignBofO,
                    const vector<int>& assignOofA,
                    const topology &Atop, 
                    bool isA, bool wang_approximation, bool use_topological_min)
{
  int N = (int)assignAofO.size();
  vector<int> Aphantoms; // A's atoms which is phantom in state B, numbered in O
  vector<int> Aphantoms_A; // numbered in A

  for(int i = 0; i < N; ++i) {
    if(assignBofO[i] == -1) {
      Aphantoms.push_back(i);
      Aphantoms_A.push_back(assignAofO[i]);
    }
  }

  // Get adjacency list to find "parent" atoms
  vector<vector<int> > adj;
  Atop.convert_bonds_to_adj_list(adj);

  // Get distances from non-phantom atoms
  vector<int> distance_A =
    calculate_int_distance_nonphantom(Atop, Aphantoms_A);

  // z-matrix like connectivity information. Note that it's indexed in "O" notation.
  vector<int> bond_atom(N, -1), angle_atom(N, -1), dihedral_atom(N, -1);

  for(int p: Aphantoms) {
    int p_of_a = assignAofO[p];

    // This conn is indexed in A
    array<int, 3> conn =
      find_zmat_like_connectivity(adj,
                                  Aphantoms_A,
                                  distance_A,
                                  Atop.types,
                                  p_of_a);
    
    int bond_o = assignOofA[conn[0]];
    int angle_o = assignOofA[conn[1]];
    int dihed_o = assignOofA[conn[2]];
    assert(bond_o >= 0);
    assert(angle_o >= 0);
    assert(dihed_o >= 0);
    bond_atom[p] = bond_o;
    angle_atom[p] = angle_o;
    dihedral_atom[p] = dihed_o;
  }

  // Here we treat atoms that exist in A and phantom in state B.
  // Thus the force constant is multiplied as in the following factors.
  double factorA = isA ? 0.0 : 1.0;
  double factorB = isA ? 1.0 : 0.0;

  if(!wang_approximation) {
    Ofs << "[ bonds ]" << endl;
    
    for(int p: Aphantoms) {
      int pa = assignAofO[p];
      int bond_o = bond_atom[p];
      int ba = assignAofO[bond_o];

      double d = (Acoords.col(pa) - Acoords.col(ba)).norm();
      if(use_topological_min) {
        d = topological_min_dist(Atop, pa, ba);
      }
      Ofs << setw(6) << p + 1 << " "
          << setw(6) << bond_o + 1 << " "
          << setw(4) << 6 << " " // Harmonic, non-excl
          << setw(8) << d * nm_of_angstrom << " "
          << setw(8) << factorA * force_bond << " "
          << setw(8) << d * nm_of_angstrom << " "
          << setw(8) << factorB * force_bond << " "
          << "; (dummy conn.) " 
          << setw(3) << Atop.resids[pa] << " "
          << setw(5) << Atop.names[pa] << " - "
          << setw(3) << Atop.resids[ba] << " "
          << setw(5) << Atop.names[ba]
          << endl;
    }
  }

  if(!wang_approximation) {
    Ofs << "[ angles ]" << endl;

    for(int p: Aphantoms) {
      int pa = assignAofO[p];
      int bond_o = bond_atom[p];
      int ba = assignAofO[bond_o];
      int angle_o = angle_atom[p];
      int aa = assignAofO[angle_o];

      double ang_rad = angle(Acoords.col(pa), Acoords.col(ba), Acoords.col(aa));
      if(use_topological_min) {
        ang_rad = topological_min_angle(Atop, pa, ba, aa) * M_PI / 180.0;
      }
      Ofs << setw(6) << p + 1 << " "
          << setw(6) << bond_o + 1 << " "
          << setw(6) << angle_o + 1 << " "
          << setw(4) << 1 << " " // harmonic
          << setw(8) << ang_rad * 180.0 / M_PI << " "
          << setw(8) << factorA * force_angluar << " "
          << setw(8) << ang_rad * 180.0 / M_PI << " "
          << setw(8) << factorB * force_angluar << " "
          << "; (dummy conn.) " 
          << setw(3) << Atop.resids[pa] << " "
          << setw(5) << Atop.names[pa] << "-"
          << setw(5) << Atop.names[ba] << "-"
          << setw(5) << Atop.names[aa]
          << endl;
    }
  }

  Ofs << "[ dihedrals ]" << endl;

  for(int p: Aphantoms) {
    int pa = assignAofO[p];
    int bond_o = bond_atom[p];
    int ba = assignAofO[bond_o];
    int angle_o = angle_atom[p];
    int aa = assignAofO[angle_o];
    int dihed_o = dihedral_atom[p];
    int da = assignAofO[dihed_o];

    // Here I assume that, there are no three atoms that are linear in state B and non-linear in state A - are there any major counterexample?
    if(unfavored_for_dihedral_restraint(Acoords, pa, ba, aa) ||
       unfavored_for_dihedral_restraint(Acoords, ba, aa, da)){
      if(verbose) {
        cerr << "Skipping dihedral restraint " << (p + 1) << "-" << (bond_o + 1) << "-" << (angle_o + 1) << "-" << (dihed_o + 1) << " because some atoms are linear" << endl;
      }
      continue;
    }

    double dih_rad = dihedral(Acoords.col(pa), Acoords.col(ba), Acoords.col(aa), Acoords.col(da));
    if(use_topological_min) {
      vector<double> dih_cands = topological_min_dihed(Atop, pa, ba, aa, da);
      assert(dih_cands.size() > 0);
      if(verbose) {
        cerr << "DEBUG Topomin dihed " << (p + 1) << " " << (bond_o + 1) << " " << (angle_o + 1) << " " << (dihed_o + 1) << " base " << dih_rad << endl;
        for(double d: dih_cands) {
          cerr << d << " ";
        }
        cerr << endl;
      }
      double dih_best_angle = dih_cands[0];
      double curbest = 1e9;
      for(double d: dih_cands) {
        double drad = d * M_PI / 180.0;
        double diff = drad - dih_rad;
        diff -= 2 * M_PI * round(diff / (2 * M_PI)); // normalize to -PI .. PI
        diff = abs(diff);
        if(curbest > diff){
          curbest = diff;
          dih_best_angle = drad;
        }
      }
      dih_rad = dih_best_angle;
    }
    Ofs << setw(6) << p + 1 << " "
        << setw(6) << bond_o + 1 << " "
        << setw(6) << angle_o + 1 << " "
        << setw(6) << dihed_o + 1 << " "
        << setw(4) << 2 << " " // "improper" to constrain angle with harmonic
        << setw(8) << dih_rad * 180.0 / M_PI << " "
        << setw(8) << factorA * force_dihedral << " "
        << setw(8) << dih_rad * 180.0 / M_PI << " "
        << setw(8) << factorB * force_dihedral << " "
        << "; (dummy conn.) " 
        << setw(3) << Atop.resids[pa] << " "
        << setw(5) << Atop.names[pa] << "-"
        << setw(5) << Atop.names[ba] << "-"
        << setw(5) << Atop.names[aa] << "-"
        << setw(5) << Atop.names[da]
        << endl;
  }
  
}

static void find_exclusions(const topology &top,
                            const vector<int> &assignOofX,
                            vector<vector<int> >& exclusions)
{
  assert(exclusions.size() != 0);
  int nexcl = top.nexcl;

  // construct bond table
  vector<vector<int> > adj;
  top.convert_bonds_to_adj_list(adj);

  for(int i = 0; i < (int)adj.size(); ++i) {
    vector<bool> visited(adj.size(), false);
    queue<pair<int, int> > bfsq;
    bfsq.emplace(i, 0);
    while(!bfsq.empty()) {
      auto e = bfsq.front();
      int v = e.first;
      int depth = e.second;
      bfsq.pop();
      if(visited[v]) continue;
      if(depth > 0) {
        exclusions[assignOofX[i]].push_back(assignOofX[v]);
      }
      visited[v] = true;
      if(depth < nexcl) {
        for(int x: adj[v]) {
          if(visited[x]) continue;
          bfsq.emplace(x, depth + 1);
        }
      }
    }
  }
}

static void make_residue_index(const pdb& pdb,
                               vector<int> *residue_indices, 
                               int *Nresidue)
{
  residue_indices->clear();

  int Nresidue_ = -1;
  pair<int, char> prev(0, ' ');
  for(size_t i = 0; i < pdb.get_numatoms(); ++i) {
    pair<int, char> newres = make_pair(pdb.get_resids()[i], pdb.get_chains()[i]);
    if(i == 0 || prev != newres) {
      prev = newres;
      Nresidue_++;
    }
    residue_indices->push_back(Nresidue_);
  }
  *Nresidue = Nresidue_;
}

static void merge_assigns(const vector<int>& assignBofA,
                          const vector<int>& assignAofB,
                          const vector<int>& residue_indices_A,
                          const vector<int>& residue_indices_B,
                          vector<int>& assignOofA,
                          vector<int>& assignOofB,
                          vector<int>& assignAofO,
                          vector<int>& assignBofO,
                          int *N)
{
  int Bunassigned = 0;
  for(int i = 0; i < (int)assignAofB.size(); ++i) {
    if(assignAofB[i] == -1) Bunassigned++;
  }
  *N = assignBofA.size() + Bunassigned;
  assignOofA = vector<int>(assignBofA.size(), -1);
  assignOofB = vector<int>(assignAofB.size(), -1);
  assignAofO = vector<int>(*N, -1);
  assignBofO = vector<int>(*N, -1);

  for(int i = 0; i < (int)assignOofA.size(); ++i) {
    assignOofA[i] = i;
    assignAofO[i] = i;
  }
  for(int j = 0, ptr = (int)assignBofA.size();
      j < (int)assignOofB.size(); ++j) {
    if(assignAofB[j] == -1) {
      assignOofB[j] = ptr;
      ptr++;
    }else{
      assignOofB[j] = assignOofA[assignAofB[j]];
    }
    assignBofO[assignOofB[j]] = j;
  }

  // Initial map was done, then reorder by the residue index 
  vector<int> oremap;
  for(int i = 0; i < *N; ++i) {
    oremap.push_back(i);
  }
  // sort by residue index
  std::stable_sort(oremap.begin(), oremap.end(), 
                   [&](const int &lhs, const int &rhs){
                      if(assignAofO[lhs] != -1 && assignAofO[rhs] != -1) {
                        int rl = residue_indices_A[assignAofO[lhs]];
                        int rr = residue_indices_A[assignAofO[rhs]];
                        if(rl != rr) { return rl < rr; }

                        // check residue B status
                        if(assignBofO[lhs] == -1) {
                          return false;
                        }else if(assignBofO[rhs] == -1) {
                          return true;
                        }else{
                          return residue_indices_B[assignBofO[lhs]] < residue_indices_B[assignBofO[rhs]];
                        }
                      }else if(assignBofO[lhs] != -1 && assignBofO[rhs] != -1) {
                        int rl = residue_indices_B[assignBofO[lhs]];
                        int rr = residue_indices_B[assignBofO[rhs]];
                        if(rl != rr) { return rl < rr; }

                        // check residue A status
                        if(assignAofO[lhs] == -1) {
                          return false;
                        }else if(assignAofO[rhs] == -1) {
                          return true;
                        }else{
                          assert(!"This should not happen as the second big else-if clause is called only when either of assignAofO == -1");
                        }
                      }else if(assignAofO[lhs] == -1) {
                        // residue index(A, B): lhs -> (-1, B) rhs -> (A, -1)
                        assert(assignBofO[rhs] == -1);
                        return false;
                      }else if(assignAofO[rhs] == -1) {
                        // lhs -> (A, -1) rhs-> (-1, B)
                        assert(assignBofO[lhs] == -1);
                        return true;
                      }else{
                        assert(!"Unexpected resindex pattern");
                      }
                   });

  vector<int> assignAofO_new, assignBofO_new;
  for(int i = 0; i < *N; ++i) {
    assignAofO_new.push_back(assignAofO[oremap[i]]);
    assignBofO_new.push_back(assignBofO[oremap[i]]);
  }
  assignAofO.swap(assignAofO_new);
  assignBofO.swap(assignBofO_new);
  vector<int> assignOofA_new(assignBofA.size(), -1);
  vector<int> assignOofB_new(assignAofB.size(), -1);
  assignOofA.swap(assignOofA_new);
  assignOofB.swap(assignOofB_new);
  for(int i = 0; i < *N; ++i) {
    if(assignAofO[i] != -1) {
      assignOofA[assignAofO[i]] = i;
    }
    if(assignBofO[i] != -1) {
      assignOofB[assignBofO[i]] = i;
    }
  }

}

static void parse_push_types_helper(const string& l,
                                    const string& types_name,
                                    int nmatch,
                                    vector<int> *funcs,
                                    vector<tuple<int, int, string> > *outputs)
{
  istringstream is(l);
  string atype;
  int Xcounts = 0;
  int func;
  for(int i = 0; i < nmatch; ++i) {
    is >> atype;
    if(atype == "X") Xcounts++;
  }
  is >> func;
  if(!is) {
    cerr << "Failed to parse [ " << types_name << " ] line: \"" << l << "\"" << endl;
    exit(1);
  }
  funcs->push_back(func);
  outputs->push_back(make_tuple(func, Xcounts, l));
}

template<typename T>
static void parse_push_types(const T& Atypes,
                             const T& Btypes,
                             const string& types_name,
                             int nmatch,
                             vector<int> *funcs,
                             vector<tuple<int, int, string> > *outputs)
{
  for(const auto& v: Atypes) {
    for(const auto &l: v.second) {
      parse_push_types_helper(l, types_name, nmatch, funcs, outputs);
    }
    for(const auto& v: Btypes) {
      if(Atypes.count(v.first) == 0) {
        for(const auto &l: v.second) {
          parse_push_types_helper(l, types_name, nmatch, funcs, outputs);
        }
      }
    }
  }
}
 
void correct_assign_by_exclusion(const topology &Atop,
                                 const topology &Btop,
                                 const pdb &Apdb,
                                 const pdb &Bpdb,
                                 vector<int> &assignBofA,
                                 vector<int> &assignAofB,
                                 vector<int> &Adepth)
{
  vector<vector<int> > Apairs_A(assignBofA.size()), Bpairs_B(assignAofB.size());
  Atop.convert_pairs_to_adj_list(Apairs_A);
  Btop.convert_pairs_to_adj_list(Bpairs_B);

  vector<int> Aresindex, Bresindex;
  int NAres, NBres;
  make_residue_index(Apdb, &Aresindex, &NAres);
  make_residue_index(Bpdb, &Bresindex, &NBres);

  bool found;
  do {
    found = false;
    int maxdepth = -1;
    int Amax, Bmax = -1;
    vector<int> assignOofA, assignOofB, assignAofO, assignBofO;
    int N;
    // try to generate current assign
    merge_assigns(assignBofA, assignAofB, 
                  Aresindex, Bresindex,
                  assignOofA, assignOofB,
                  assignAofO, assignBofO,
                  &N);
    // generate exclusion
    vector<vector<int> > Aexclusions(N), Bexclusions(N);
    find_exclusions(Atop, assignOofA, Aexclusions);
    find_exclusions(Btop, assignOofB, Bexclusions);

    // generate pairs and convert
    vector<vector<int> > Apairs(N), Bpairs(N);
    for(int i = 0; i < (int) Apairs_A.size(); ++i) {
      for(int v: Apairs_A[i]) {
        Apairs[assignOofA[i]].push_back(assignOofA[v]);
      }
    }
    for(int i = 0; i < (int) Bpairs_B.size(); ++i) {
      for(int v: Bpairs_B[i]) {
        Bpairs[assignOofB[i]].push_back(assignOofB[v]);
      }
    }
    
    // If a pair exists in both states, it's OK.
    // If a pair exists in neither states, it's OK.
    // Otherwise, such a pair should not appear in the final assignment
    vector<int> xors;
    for(int i = 0; i < (int)assignOofA.size(); ++i) {
      if(assignAofO[i] == -1 ||
         assignBofO[i] == -1){
        continue;
      }
      for(int mode = 0; mode < 2; ++mode) {
        xors.clear();
        if(mode == 0) {
          // excls
          sort(Aexclusions[i].begin(), Aexclusions[i].end());
          sort(Bexclusions[i].begin(), Bexclusions[i].end());
          set_symmetric_difference(Aexclusions[i].begin(),
                                   Aexclusions[i].end(),
                                   Bexclusions[i].begin(),
                                   Bexclusions[i].end(),
                                   back_inserter(xors));
        }else if(mode == 1) {
          sort(Apairs[i].begin(), Apairs[i].end());
          sort(Bpairs[i].begin(), Bpairs[i].end());
          set_symmetric_difference(Apairs[i].begin(),
                                   Apairs[i].end(),
                                   Bpairs[i].begin(),
                                   Bpairs[i].end(),
                                   back_inserter(xors));
        }else{
          assert(false);
        }
        for(int j: xors) {
          if(assignAofO[j] == -1 ||
             assignBofO[j] == -1){
            continue;
          }
          // found unmatching exclusion
          found = true;
          int d;
          d = Adepth[assignAofO[i]];
          if(d > maxdepth) {
            maxdepth = d;
            Amax = assignAofO[i];
            Bmax = assignBofO[i];
          }
          d = Adepth[assignAofO[j]];
          if(d > maxdepth) {
            maxdepth = d;
            Amax = assignAofO[j];
            Bmax = assignBofO[j];
          }
        }
      }
    }
    if(found) {
      cerr << "Found unmatching exclusion" << endl;
      cerr << "Removing atom " 
           << Atop.resnames[Amax] << "_" 
           << Atop.resids[Amax] << ":"
           << Atop.names[Amax] << " - "
           << Btop.resnames[Bmax] << "_" 
           << Btop.resids[Bmax] << ":"
           << Btop.names[Bmax]
           << " from the list of matched atoms" << endl;
      assignBofA[Amax] = -1;
      assignAofB[Bmax] = -1;
    }
  }while(found);
}

// linear interpolation
// x = 0 -> a
// x = 1 -> b
double li(double x, double a, double b)
{
  return a * (1 - x) + b * x;
}


static void do_check_cmap_add_zerofill(const topology &Atop, const topology &Btop,
                                       const vector<int> &assignOofA, const vector<int> &assignBofO,
                                       map<topology::cmaptype, tuple<int, int, vector<double>>> *cmaptypes_plus0)
{
  assert(cmaptypes_plus0->size() > 0);
  for(const auto& cmA: Atop.cmaps) {
    int a = get<0>(cmA);
    int b = get<1>(cmA);
    int c = get<2>(cmA);
    int d = get<3>(cmA);
    int e = get<4>(cmA);
    int func = get<5>(cmA);

    string at = Atop.bondatomtypes.at(Atop.types[a]);
    string bt = Atop.bondatomtypes.at(Atop.types[b]);
    string ct = Atop.bondatomtypes.at(Atop.types[c]);
    string dt = Atop.bondatomtypes.at(Atop.types[d]);
    string et = Atop.bondatomtypes.at(Atop.types[e]);

    assert(assignOofA[a] != -1);
    assert(assignOofA[b] != -1);
    assert(assignOofA[c] != -1);
    assert(assignOofA[d] != -1);
    assert(assignOofA[e] != -1);

    int a_in_B = assignBofO[assignOofA[a]];
    int b_in_B = assignBofO[assignOofA[b]];
    int c_in_B = assignBofO[assignOofA[c]];
    int d_in_B = assignBofO[assignOofA[d]];
    int e_in_B = assignBofO[assignOofA[e]];

    string at_in_B = (a_in_B != -1 ? Btop.bondatomtypes.at(Btop.types[a_in_B]) : "PHA");
    string bt_in_B = (b_in_B != -1 ? Btop.bondatomtypes.at(Btop.types[b_in_B]) : "PHA");
    string ct_in_B = (c_in_B != -1 ? Btop.bondatomtypes.at(Btop.types[c_in_B]) : "PHA");
    string dt_in_B = (d_in_B != -1 ? Btop.bondatomtypes.at(Btop.types[d_in_B]) : "PHA");
    string et_in_B = (e_in_B != -1 ? Btop.bondatomtypes.at(Btop.types[e_in_B]) : "PHA");

    // TODO FIXME: as a pathological case, number of discretization may not match throughout entries.

    topology::cmaptype keyA = std::make_tuple(at, bt, ct, dt, et, func);
    topology::cmaptype keyB = std::make_tuple(at_in_B, bt_in_B, ct_in_B, dt_in_B, et_in_B, func);
    if(cmaptypes_plus0->find(keyB) ==
       cmaptypes_plus0->end()) {
      if(verbose) {
        cout << "Adding 0-filled CMAP entry for "
             << at_in_B << "-"
             << bt_in_B << "-"
             << ct_in_B << "-"
             << dt_in_B << "-"
             << et_in_B << " (func " << func;
      }
      const auto& baseentry = Atop.cmaptypes.at(keyA);
      int phi_split = get<0>(baseentry);
      int psi_split = get<1>(baseentry);
      vector<double> newentry(get<2>(baseentry).size(), 0.);
      if(verbose) {
        cout << ", numentry = " << get<2>(baseentry).size() << ")" << endl;
      }
      (*cmaptypes_plus0)[keyB] = std::make_tuple(phi_split, psi_split, newentry);
    }
  }
}

// Hey, could you clean up this mess?
void generate_topology(const string& outfilename,
                       const vector<string>& argv,
                       const string& title,
                       const topology &Atop,
                       const topology &Btop,
                       const Matrix3Xd &Acoords,
                       const Matrix3Xd &Bcoords,
                       const vector<double>& stateA_weights,
                       const vector<double>& stateB_weights,
                       const vector<int> &assignOofA,
                       const vector<int> &assignOofB, 
                       const vector<int> &assignAofO,
                       const vector<int> &assignBofO,
                       const struct dummy_anchoring_parameters& cparams,
                       double rest_weighting,
                       const vector<string> &rest_exclude_atomtypes,
                       bool gen_exclusion,
                       bool wang_approximation,
                       bool disable_pairs_error,
                       bool use_topological_min,
                       const string& restraints_atomname)
{
  int N = assignAofO.size();

  // output 
  ofstream Ofs(outfilename.c_str());
  Ofs.setf(ios::scientific);
  Ofs.setf(ios::right);
  Ofs << setprecision(5);

  Ofs << "; generated by:" << endl;
  Ofs << ";  " << argv[0];
  for(int i = 1; i < (int)argv.size(); ++i) {
    Ofs << " " << argv[i];
  }
  Ofs << endl;
  
  // Defaults section
  if(Atop.defaults != Btop.defaults ||
     !(Atop.defaults == topology::topotype::AMBER ||
       Atop.defaults == topology::topotype::CHARMM ||
       Atop.defaults == topology::topotype::OPLS)) {
    cerr << "Unsupported defaults type" << endl;
    exit(1);
  }
  Ofs << "[ defaults ]" << endl;
  Ofs << "; nbfunc        comb-rule       gen-pairs       fudgeLJ fudgeQQ" << endl;
  if(Atop.defaults == topology::topotype::AMBER) {
    Ofs << "1               2               yes             0.5     0.8333" << endl;
  }else if(Atop.defaults == topology::topotype::CHARMM) {
    Ofs << "1 2 yes 1.0 1.0" << endl;
  }else if(Atop.defaults == topology::topotype::OPLS) {
    Ofs << "1 3 yes 0.5 0.5" << endl;
  }else{
    throw std::runtime_error("Unsupported defaults type");
  }
  Ofs << endl;

  // atomtypes section
  vector<string> atomtypes;
  Ofs << "[ atomtypes ]" << endl;
  for(const auto& v: Atop.atomtypes) {
    Ofs << v.second << endl;
    {
      istringstream is(v.second);
      string at;
      is >> at;
      atomtypes.push_back(at);
    }
  }
  for(const auto& v: Btop.atomtypes) {
    if(Atop.atomtypes.count(v.first) == 0) {
      Ofs << v.second << endl;
      {
        istringstream is(v.second);
        string at;
        is >> at;
        atomtypes.push_back(at);
      }
    }
  }
  // phantom type
  Ofs << "PHA PHA 0.0000 0.0000 A 0.0000 0.0000" << endl;
  Ofs << endl;
  atomtypes.push_back("PHA");

  // pairtypes section
  if(Atop.pairtypes.size() + 
     Btop.pairtypes.size() > 0) {
    Ofs << "[ pairtypes ]" << endl;
    
    for(const auto& v: Atop.pairtypes) {
      Ofs << v.second << endl;
    }
    for(const auto& v: Btop.pairtypes) {
      if(Atop.pairtypes.count(v.first) == 0) {
        Ofs << v.second << endl;
      }
    }
    Ofs << endl;
  }

  // bondtypes
  if(Atop.bondtypes.size() +
     Btop.bondtypes.size() > 0) {
    Ofs << "[ bondtypes ]" << endl;

    vector<int> funcs;
    vector<tuple<int, int, string> > outputs;

    parse_push_types(Atop.bondtypes,
                     Btop.bondtypes,
                     "bondtypes",
                     2,
                     &funcs,
                     &outputs);
   
    sort(funcs.begin(), funcs.end());
    funcs.erase(unique(funcs.begin(), funcs.end()), funcs.end());
    
    for(int f: funcs) {
      // FIXME: does not suppor morse potential
      for(const string& at: atomtypes) {
        ostringstream os;
        os << "PHA " << at << " " << f << " 0.0 0.0";
        outputs.push_back(make_tuple(f, 99, os.str())); // wanna put at the end
      }
    }
    sort(outputs.begin(), outputs.end());
    for(const auto& tps: outputs) {
      Ofs << get<2>(tps) << endl;
    }
    Ofs << endl;
  }

  // angletypes
  if(Atop.angletypes.size() +
     Btop.angletypes.size() > 0) {
    Ofs << "[ angletypes ]" << endl;

    vector<int> funcs;
    vector<tuple<int, int, string> > outputs;

    parse_push_types(Atop.angletypes,
                     Btop.angletypes,
                     "angletypes",
                     3,
                     &funcs,
                     &outputs);

    sort(funcs.begin(), funcs.end());
    funcs.erase(unique(funcs.begin(), funcs.end()), funcs.end());
    
    for(int f: funcs) {
      int zeros = 2;
      if(f == 5) { zeros = 4; }
      
      for(const string& at: atomtypes) {
        for(const string& at2: atomtypes) {
          ostringstream os;
          os << "PHA " << at << " " << at2 << " " << f;
          for(int i = 0; i < zeros; ++i) {
            os << " 0.0";
          }
          outputs.push_back(make_tuple(f, 99, os.str()));
        }
      }

      for(const string& at: atomtypes) {
        if(at == "PHA") continue;
        for(const string& at2: atomtypes) {
          ostringstream os;
          os << at << " PHA " << at2 << " " << f;
          for(int i = 0; i < zeros; ++i) {
            os << " 0.0";
          }
          outputs.push_back(make_tuple(f, 99, os.str()));
        }
      }
    }
    sort(outputs.begin(), outputs.end());
    for(const auto& tps: outputs) {
      Ofs << get<2>(tps) << endl;
    }

    Ofs << endl;
  }

  // dihedraltypes
  // to combat gromacs issue #1901, and to match proper dihedrals,
  // first list non-improper dihedrals, 
  // then list improper diherals.
  // Note that multiplicity matters. 
  // [ dihedraltypes ]
  // A B C D 9 ... 1
  // A B C D 9 ... 2
  // D E F G 9 ...
  // X A B X 9 ...
  // A B D E 2 ...
  // X X A D 2 ...
  if(Atop.dihedraltypes.size() +
     Btop.dihedraltypes.size() > 0) {
    Ofs << "[ dihedraltypes ]" << endl;

    // func, multiplicity
    vector<pair<int, int> > funcs;
    // i do know this is a dirrrrrrrrty hack
    // (is_2, has_X) -> inputs
    map<pair<bool, bool>, vector<string> > outputs;
    
    for(const auto& v: Atop.dihedraltypes) {
      for(const auto &l: v.second) {
        istringstream is(l);
        string a, b, c, d;
        int func;
        int multi = -1;
        is >> a >> b >> c >> d >> func;
        if(func == 1 || func == 4 || func == 9) {
          double dummy;
          is >> dummy >> dummy >> multi;
        }
        if(!is) {
          cerr << "Failed to parse dihedraltypes line: \"" << l << "\"" << endl;
          exit(1);
        }
        funcs.push_back(make_pair(func, multi));
        bool has_x = (a == "X" || b == "X" || c == "X" || d == "X");
        bool is_2 = func == 2;
        outputs[make_pair(is_2, has_x)].push_back(l);
      }
    }
    for(const auto& v: Btop.dihedraltypes) {
      if(Atop.dihedraltypes.count(v.first) == 0) {
        for(const auto &l: v.second) {
          istringstream is(l);
          string a, b, c, d;
          int func;
          int multi = -1;
          is >> a >> b >> c >> d >> func;
          if(func == 1 || func == 4 || func == 9) {
            double dummy;
            is >> dummy >> dummy >> multi;
          }
          if(!is) {
            cerr << "Failed to parse dihedraltypes line: \"" << l << "\"" << endl;
            exit(1);
          }
          funcs.push_back(make_pair(func, multi));
          bool has_x = (a == "X" || b == "X" || c == "X" || d == "X");
          bool is_2 = func == 2;
          outputs[make_pair(is_2, has_x)].push_back(l);
        }
      }
    }
    
    sort(funcs.begin(), funcs.end());
    funcs.erase(unique(funcs.begin(), funcs.end()), funcs.end());
    
    for(auto p: funcs) {
      int f = p.first;
      int m = p.second;
      int zeros = 2;
      if(f == 3) { zeros = 6; }
      if(f == 5) { zeros = 4; }
      
      {
        ostringstream os;
        os << "PHA X X X " << f;
        for(int i = 0; i < zeros; ++i) {
          os << " 0.0";
        }
        if(m >= 0) {
          os << " " << m;
        }
        outputs[make_pair(f, true)].push_back(os.str());
      }
    }

    for(auto p: funcs) {
      int f = p.first;
      int m = p.second;
      int zeros = 2;
      if(f == 3) { zeros = 6; }
      if(f == 5) { zeros = 4; }

      {
        ostringstream os;
        os << "X PHA X X " << f;
        for(int i = 0; i < zeros; ++i) {
          os << " 0.0";
        }
        if(m >= 0) {
          os << " " << m;
        }
        outputs[make_pair(f, true)].push_back(os.str());
      }
    }

    outputs[make_pair(2, true)].push_back("X X PHA X 2 0.0 0.0");
    outputs[make_pair(2, true)].push_back("X X X PHA 2 0.0 0.0");

    for(const auto& p: outputs){
      for(const auto& l: p.second) {
        Ofs << l << endl;
      }
    }
    Ofs << endl;
  }

  // cmaptypes.
  // During the FEP, one of state A and B may not have the cmap entry for the perturbed bond-atom type. In such a case, we need to supply an entry of the cmaptypes with all values being 0.
  if(Atop.cmaptypes.size() > 0 || Btop.cmaptypes.size() > 0) {
    map<topology::cmaptype, tuple<int, int, vector<double>>> cmaptypes_plus0(Atop.cmaptypes);
    // We assume cmaptypes are identical between two topologies.
    // This should be valid as long as we use like #include "charmmxxff/cmap.itp" ...
    if(Atop.cmaptypes.size() != Btop.cmaptypes.size()) {
      throw runtime_error("[ cmaptypes ] entries are not identical between two states (size check failed)");
    }
    for(const auto& ceA: Atop.cmaptypes) {
      map<topology::cmaptype, tuple<int, int, vector<double>>>::const_iterator it;
      // in fact we can just use std::equal and lambda, but to make more informative error message...
      if((it = Btop.cmaptypes.find(ceA.first)) != Btop.cmaptypes.end()) {
        // entry exists in Btop
        bool identical_entry =
          (get<0>(ceA.second) == get<0>(it->second)) &&
          (get<1>(ceA.second) == get<1>(it->second)) &&
          (get<2>(ceA.second).size() == get<2>(it->second).size()) &&
          std::equal(get<2>(ceA.second).begin(), get<2>(ceA.second).end(), get<2>(it->second).begin());
        if(!identical_entry) {
          cerr << "CMAP error: "
               << get<0>(ceA.first) << "-"
               << get<1>(ceA.first) << "-"
               << get<2>(ceA.first) << "-"
               << get<3>(ceA.first) << "-"
               << get<4>(ceA.first) << " (functype "
               << get<5>(ceA.first) << ")" << endl;
          throw runtime_error("[ cmaptypes ] entries are not identical between two states");
        }
      }else{
        cerr << "CMAP error: "
               << get<0>(ceA.first) << "-"
               << get<1>(ceA.first) << "-"
               << get<2>(ceA.first) << "-"
               << get<3>(ceA.first) << "-"
               << get<4>(ceA.first) << " (functype "
               << get<5>(ceA.first) << ")" << endl;
        throw runtime_error("[ cmaptypes ] entry in state A does not exist in state B");
      }
    }
    // With the equal size, iterating over A should be sufficient.

    // Next: check CMAP usage, and add 0-filled entry if necessray.
    do_check_cmap_add_zerofill(Atop, Btop, assignOofA, assignBofO, &cmaptypes_plus0);
    do_check_cmap_add_zerofill(Btop, Atop, assignOofB, assignAofO, &cmaptypes_plus0);

    Ofs << "[ cmaptypes ]" << endl;
    for(const auto &entry: cmaptypes_plus0) {
      const topology::cmaptype &key = entry.first;
      const auto &value = entry.second;
      int phi_split = get<0>(value);
      int psi_split = get<1>(value);
      const vector<double>& potentials = get<2>(value);
      Ofs << get<0>(key) << " "
          << get<1>(key) << " "
          << get<2>(key) << " "
          << get<3>(key) << " "
          << get<4>(key) << " "
          << get<5>(key) << " "
          << phi_split << " "
          << psi_split << "\\"
          << endl;
      Ofs << std::fixed;
      Ofs.precision(8);
      for(int i = 0; i < phi_split; ++i) {
        for(int j = 0; j < psi_split; ++j) {
          Ofs << potentials[i * psi_split + j] << (j == psi_split - 1 ? "" : " ");
        }
        Ofs << (i == phi_split - 1 ? "" : "\\") << endl;
      }
    }
  }

  if(!Atop.nbfixes.empty()) {
    Ofs << "[ nonbond_params ]" << endl;
    
    bool identical =
      Atop.nbfixes.size() == Btop.nbfixes.size() &&
      std::equal(Atop.nbfixes.begin(), Atop.nbfixes.end(), Btop.nbfixes.begin());
    if(!identical) {
      cerr << "[ nonbond_params ] section are not identical between two topologies" << endl;
      throw std::runtime_error("Reimplement saner nonbond_params process to topology.cpp");
    }
    for(const string& l: Atop.nbfixes) {
      Ofs << l << endl;
    }
  }

  // moleculetype section
  Ofs << "[ moleculetype ]" << endl;
  Ofs << ";name  nrexcl" << endl;
  {
    int nexcl = Atop.nexcl;
    if(gen_exclusion) {
      nexcl = 0;
    }
    Ofs << "merged " << nexcl << endl;
  }
  Ofs << endl;

  // atoms section
  Ofs << "[ atoms ]" << endl;
  Ofs << ";   nr  type  resi  res  atom  cgnr     charge      mass       typeB chargeB massB" << endl;

  for(int i = 0; i < N; ++i) {
    double mass;
    if(assignAofO[i] != -1 && assignBofO[i] != -1) {
      mass = std::max(Atop.masses[assignAofO[i]], Btop.masses[assignBofO[i]]);
    }else if(assignAofO[i] != -1) {
      mass = Atop.masses[assignAofO[i]];
    }else{
      mass = Btop.masses[assignBofO[i]];
    }
    
    // Fill in A-state
    if(assignAofO[i] != -1) {
      int aptr = assignAofO[i];
      Ofs << setw(4) << (i + 1)
          << " " << setw(3) << Atop.types[aptr]
          << " " << setw(4) << Atop.resids[aptr]
          << " " << setw(5) << Atop.resnames[aptr]
          << " " << setw(6) << Atop.names[aptr]
          << " " << setw(4) << (i + 1) // cgnr
          << " " << setw(12) << Atop.charges[aptr]
          << " " << setw(12) << mass;
    }else{
      int bptr = assignBofO[i];
      assert(bptr != -1);
      Ofs << setw(4) << (i + 1)
          << " " << setw(3) << "PHA"
          << " " << setw(4) << Btop.resids[bptr]
          << " " << setw(5) << Btop.resnames[bptr]
          << " " << setw(6) << Btop.names[bptr]
          << " " << setw(4) << (i + 1) // cgnr
          << " " << setw(12) << 0.000
          << " " << setw(12) << mass;
    }

    // Fill in B-state
    if(assignBofO[i] != -1) {
      int bptr = assignBofO[i];
      Ofs << " " << setw(3) << Btop.types[bptr]
          << " " << setw(12) << Btop.charges[bptr]
          << " " << setw(12) << mass;
    }else{
      Ofs << " " << "PHA"
          << " " << setw(12) << 0.000
          << " " << setw(12) << mass;
    }
    Ofs << endl;
  }
  Ofs << endl;

  // bonds section
  Ofs << "[ bonds ]" << endl;
  for(const auto& v: Atop.bonds) {
    topology::bondkeytype k = v.first;
    int x = std::get<0>(k);
    int y = std::get<1>(k);
    int func = std::get<2>(k);
    
    int x_in_O = assignOofA[x];
    int y_in_O = assignOofA[y];
    topology::bondkeytype keyB = 
      make_tuple(assignBofO[x_in_O],
                 assignBofO[y_in_O],
                 func);
    if(func != 1 && func != 6) {
      throw runtime_error("Bond func != 1, 6 not supported");
    }
    const vector<double>& Afactors = v.second;
    vector<double> Bfactors(2, 0);
    if(Btop.bonds.count(keyB) > 0) {
      Bfactors = Btop.bonds.at(keyB);
    }else if(Afactors.size() > 0){
      Bfactors[0] = Afactors[0];
      if(wang_approximation){
        Bfactors[1] = Afactors[1];
      }else{
        Bfactors[1] = 0.0;
      }
    }else{
      Bfactors.clear();
    }
    // output A-listed bonds first
    Ofs << x_in_O + 1 << " " << y_in_O + 1 << " " 
        << func;
    for(auto v: Afactors) {
      Ofs << " " << v;
    }
    for(auto v: Bfactors) {
      Ofs << " " << v;
    }
    Ofs << endl;
  }
  // Remaining: not exists in A but exists in B
  for(const auto& v: Btop.bonds) {
    topology::bondkeytype k = v.first;
    int x = std::get<0>(k);
    int y = std::get<1>(k);
    int func = std::get<2>(k);
    
    int x_in_O = assignOofB[x];
    int y_in_O = assignOofB[y];
    topology::bondkeytype keyA = 
      make_tuple(assignAofO[x_in_O],
                 assignAofO[y_in_O],
                 func);
    if(Atop.bonds.count(keyA) > 0) {
      continue;
    }
    const vector<double>& Bfactors = v.second;
    Ofs << x_in_O + 1 << " " << y_in_O + 1 << " " 
        << func;
    if(Bfactors.size() > 0) {
      Ofs << " " << Bfactors[0]
          << " " << (wang_approximation ? Bfactors[1] : 0.00);
    }
    for(auto v: Bfactors) {
      Ofs << " " << v;
    }
    Ofs << endl;
  }
  Ofs << endl;

  // angles section (FIXME: kopipe)
  Ofs << "[ angles ]" << endl;
  for(const auto& v: Atop.angles) {
    topology::anglekeytype k = v.first;
    int x = std::get<0>(k);
    int y = std::get<1>(k);
    int z = std::get<2>(k);
    int func = std::get<3>(k);
    
    int x_in_O = assignOofA[x];
    int y_in_O = assignOofA[y];
    int z_in_O = assignOofA[z];
    topology::anglekeytype keyB = 
      make_tuple(assignBofO[x_in_O],
                 assignBofO[y_in_O],
                 assignBofO[z_in_O],
                 func);
    if(func != 1 && func != 5) {
      throw runtime_error("Angle func != {1, 5} not supported");
    }
    const vector<double>& Afactors = v.second;
    vector<double> Bfactors;
    if(Btop.angles.count(keyB) > 0) {
      Bfactors = Btop.angles.at(keyB);
    }else if(Afactors.size() > 0){
      if(wang_approximation) {
        Bfactors = Afactors;
      }else{
        if(func == 1) {
          Bfactors.resize(2);
          Bfactors[0] = Afactors[0];
          Bfactors[1] = 0.0;
        }else if(func == 5) {
          // only spring constants affected
          Bfactors.resize(4);
          Bfactors[0] = Afactors[0];
          Bfactors[1] = 0.0;
          Bfactors[2] = Afactors[2];
          Bfactors[3] = 0.0;
        }
      }
    }else{
      Bfactors.clear();
    }
    // output A-listed angles first
    Ofs << x_in_O + 1 << " " << y_in_O + 1 << " " 
        << z_in_O + 1 << " "
        << func;
    for(auto v: Afactors) {
      Ofs << " " << setw(12) << v;
    }
    for(auto v: Bfactors) {
      Ofs << " " << setw(12) << v;
    }
    Ofs << endl;
  }
  // Remaining: not exists in A but exists in B
  for(const auto& v: Btop.angles) {
    topology::anglekeytype k = v.first;
    int x = std::get<0>(k);
    int y = std::get<1>(k);
    int z = std::get<2>(k);
    int func = std::get<3>(k);
    
    int x_in_O = assignOofB[x];
    int y_in_O = assignOofB[y];
    int z_in_O = assignOofB[z];
    topology::anglekeytype keyA = 
      make_tuple(assignAofO[x_in_O],
                 assignAofO[y_in_O],
                 assignAofO[z_in_O],
                 func);
    if(Atop.angles.count(keyA) > 0) {
      continue;
    }
    const vector<double>& Bfactors = v.second;
    Ofs << x_in_O + 1 << " " << y_in_O + 1 << " " 
        << z_in_O + 1 << " " 
        << func;
    if(Bfactors.size() > 0) {
      if(func == 1) {
        Ofs << " " << setw(12) << Bfactors[0]
            << " " << setw(12) << (wang_approximation ? Bfactors[1] : 0.00);
      }else{
        Ofs << " " << setw(12) << Bfactors[0]
            << " " << setw(12) << (wang_approximation ? Bfactors[1] : 0.00)
            << " " << setw(12) << Bfactors[2]
            << " " << setw(12) << (wang_approximation ? Bfactors[3] : 0.00);
      }
    }
    for(auto v: Bfactors) {
      Ofs << " " << setw(12) << v;
    }
    Ofs << endl;
  }
  Ofs << endl;

  // dihedrals (FIXME: kopipe again)
  Ofs << "[ dihedrals ]" << endl;
  for(const auto& v: Atop.diheds) {
    topology::dihedkeytype k = v.first;
    int x = std::get<0>(k);
    int y = std::get<1>(k);
    int z = std::get<2>(k);
    int w = std::get<3>(k);
    int func = std::get<4>(k);
    int addendum = std::get<5>(k);
    
    int x_in_O = assignOofA[x];
    int y_in_O = assignOofA[y];
    int z_in_O = assignOofA[z];
    int w_in_O = assignOofA[w];
    topology::dihedkeytype keyB = 
      make_tuple(assignBofO[x_in_O],
                 assignBofO[y_in_O],
                 assignBofO[z_in_O],
                 assignBofO[w_in_O],
                 func,
                 addendum);
    if(!(func == 1 || func == 2 || func == 3 || func == 4 || func == 9)) {
      throw runtime_error("dihed func not in {1, 2, 3, 4, 9} not supported");
    }
    const vector<double>& Afactors = v.second;
    vector<double> Bfactors(Afactors);
    if(Btop.diheds.count(keyB) > 0) {
      Bfactors = Btop.diheds.at(keyB);
    }else{
      for(auto &&v: Bfactors) {
        v = 0.0;
      }
      if(func == 1 || func == 4 || func == 9) {
        Bfactors[0] = Afactors[0];
        Bfactors[2] = Afactors[2];
      }
    }
    // output A-listed angles first
    Ofs << x_in_O + 1 << " " << y_in_O + 1 << " " 
        << z_in_O + 1 << " " << w_in_O + 1 << " "
        << func;
    for(size_t i = 0; i < Afactors.size(); ++i) {
      double v = Afactors[i];
      if((func == 1 || func == 4 || func == 9) && i == 2) {
        Ofs << " " << setw(2) << (int) v;
      }else{
        Ofs << " " << setw(12) << v;
      }
    }
    for(size_t i = 0; i < Bfactors.size(); ++i) {
      double v = Bfactors[i];
      if((func == 1 || func == 4 || func == 9) && i == 2) {
        Ofs << " " << setw(2) << (int) v;
      }else{
        Ofs << " " << setw(12) << v;
      }
    }
    Ofs << endl;
  }
  // Remaining: not exists in A but exists in B
  for(const auto& v: Btop.diheds) {
    topology::dihedkeytype k = v.first;
    int x = std::get<0>(k);
    int y = std::get<1>(k);
    int z = std::get<2>(k);
    int w = std::get<3>(k);
    int func = std::get<4>(k);
    int addendum = std::get<5>(k);
    
    int x_in_O = assignOofB[x];
    int y_in_O = assignOofB[y];
    int z_in_O = assignOofB[z];
    int w_in_O = assignOofB[w];
    topology::dihedkeytype keyA = 
      make_tuple(assignAofO[x_in_O],
                 assignAofO[y_in_O],
                 assignAofO[z_in_O],
                 assignAofO[w_in_O],
                 func,
                 addendum);
    if(Atop.diheds.count(keyA) > 0) {
      continue;
    }
    const vector<double>& Bfactors = v.second;
    Ofs << x_in_O + 1 << " " << y_in_O + 1 << " " 
        << z_in_O + 1 << " " << w_in_O + 1 << " "
        << func;
    if(func == 1 || func == 4 || func == 9) {
      Ofs << " " << setw(12) << Bfactors[0]
          << " " << setw(12) << 0.00
          << " " << setw(2) << (int)Bfactors[2];
    }else{
      for(double v: Bfactors) {
        (void) v;
        Ofs << " " << setw(12) << 0.00;
      }
    }
    for(size_t i = 0; i < Bfactors.size(); ++i) {
      double v = Bfactors[i];
      if((func == 1 || func == 4 || func == 9) && i == 2) {
        Ofs << " " << setw(2) << (int) v;
      }else{
        Ofs << " " << setw(12) << v;
      }
    }
    Ofs << endl;
  }
  Ofs << endl;

  if(!Atop.cmaps.empty() || !Atop.cmaps.empty()) {
    Ofs << "[ cmap ]" << endl;
    // for cmap, we only have to merge everything
    set<topology::cmapkeytype> merged;
    for(const auto& v: Atop.cmaps) {
      int a = std::get<0>(v);
      int b = std::get<1>(v);
      int c = std::get<2>(v);
      int d = std::get<3>(v);
      int e = std::get<4>(v);
      int func = std::get<5>(v);

      int a_in_O = assignOofA[a];
      int b_in_O = assignOofA[b];
      int c_in_O = assignOofA[c];
      int d_in_O = assignOofA[d];
      int e_in_O = assignOofA[e];

      merged.insert(make_tuple(a_in_O, b_in_O, c_in_O, d_in_O, e_in_O, func));
    }
    for(const auto& v: Btop.cmaps) {
      int a = std::get<0>(v);
      int b = std::get<1>(v);
      int c = std::get<2>(v);
      int d = std::get<3>(v);
      int e = std::get<4>(v);
      int func = std::get<5>(v);

      int a_in_O = assignOofB[a];
      int b_in_O = assignOofB[b];
      int c_in_O = assignOofB[c];
      int d_in_O = assignOofB[d];
      int e_in_O = assignOofB[e];

      merged.insert(make_tuple(a_in_O, b_in_O, c_in_O, d_in_O, e_in_O, func));
    }

    for(const auto& v: merged) {
      int a = std::get<0>(v);
      int b = std::get<1>(v);
      int c = std::get<2>(v);
      int d = std::get<3>(v);
      int e = std::get<4>(v);
      int func = std::get<5>(v);

      Ofs << a + 1 << " "
          << b + 1 << " "
          << c + 1 << " "
          << d + 1 << " "
          << e + 1 << " " << func << endl;
    }
    Ofs << endl;
  }

  // pairs
  Ofs << "[ pairs ]" << endl;
  for(int state = 0; state < 2; ++state){
    const topology &Ptop = (state == 0 ? Atop : Btop);
    const topology &Qtop = (state == 0 ? Btop : Atop);
    const vector<int> assignOofP = (state == 0 ? assignOofA : assignOofB);
    const vector<int> assignOofQ = (state == 0 ? assignOofB : assignOofA);
    const vector<int> assignPofO = (state == 0 ? assignAofO : assignBofO);
    const vector<int> assignQofO = (state == 0 ? assignBofO : assignAofO);
    
    for(const auto& v: Ptop.pairs) {
      auto k = v.first;
      int x = k.first;
      int y = k.second;
      int func = v.second;
      
      int x_in_O = assignOofP[x];
      int y_in_O = assignOofP[y];
      if(func != 1) {
        throw runtime_error("pairs func != 1 not supported");
      }
      bool pairs_perturbed = false;
      int x_in_Q = assignQofO[x_in_O];
      int y_in_Q = assignQofO[y_in_O];
      if(x_in_Q != -1 && y_in_Q != -1) {
        // x and y exist in state Q, check Q pairs
        if(Qtop.pairs.count({x_in_Q, y_in_Q}) == 0 && Qtop.pairs.count({y_in_Q, x_in_Q}) == 0) {
          // special case: 1-4 in state P turned to 1-2 or 1-3 (1-5 should not appear if exclusion check works...) in state Q
          pairs_perturbed = true;
        }else{
          if(state == 1) {
            // This case is already output in state 0
            continue;
          }
        }
      }
      if(pairs_perturbed){
        if(!disable_pairs_error){
          throw runtime_error("The case with pairs are perturbed cannot be run on unpatched GROMACS");
        }
        // Here I assume at least one of pairs are 1-4, but need sanity checks
        if(!(Ptop.nexcl >= 3 && Qtop.nexcl >= 3)) {
          throw runtime_error("nexcl should be greater than or equal to 3");
        }
        func = 3;
        
        Ofs << x_in_O + 1 << " " << y_in_O + 1 << " "
            << func << " "
            << (state == 0 ? 1.0 : 0.0) << " "
            << (state == 0 ? 0.0 : 1.0) << endl;
      }else{
        Ofs << x_in_O + 1 << " " << y_in_O + 1 << " " 
            << func << endl;
      }
    }
  }
  Ofs << endl;

  // bonded forces for phantom atoms
  {
    Ofs << "; Atoms existing in A state" << endl;
    output_dummies(Ofs, Acoords, 
                   cparams.force_bond,
                   cparams.force_angular,
                   cparams.force_dihedral,
                   assignAofO,
                   assignBofO,
                   assignOofA,
                   Atop, true, wang_approximation, use_topological_min);
    Ofs << "; Atoms existing in B state" << endl;
    output_dummies(Ofs, Bcoords, 
                   cparams.force_bond,
                   cparams.force_angular,
                   cparams.force_dihedral,
                   assignBofO,
                   assignAofO,
                   assignOofB,
                   Btop, false, wang_approximation, use_topological_min);
  }
  Ofs << endl;

  if(gen_exclusion) {
    Ofs << "[ exclusions ]" << endl;
    vector<vector<int> > exclusions(N);
    // generate intrinsic exclusions
    find_exclusions(Atop, assignOofA, exclusions);
    find_exclusions(Btop, assignOofB, exclusions);

    // generate dummy-dummy exclusions
    for(int i = 0; i < N; ++i) {
      if(assignAofO[i] == -1 ||
         assignBofO[i] == -1) {
        for(int j = 0; j < N; ++j) {
          if((assignAofO[i] == -1 &&
              assignBofO[j] == -1) ||
             (assignBofO[i] == -1 &&
              assignAofO[j] == -1)) {
            // A phantom vs B phantom
            exclusions[i].emplace_back(j);
          }
        }
      }
    }
    
    for(int i = 0; i < N; ++i) {
      vector<int> &v = exclusions[i];
      sort(v.begin(), v.end());
      v.erase(unique(v.begin(), v.end()), v.end());
      ostringstream os;
      Ofs << (i + 1);
      for(int e: v) {
        assert(e != i);
        //int iresid = -1, eresid = -1;
        string iname, ename;
        if(assignAofO[i] == -1 &&
           assignBofO[e] == -1) {
          int iid = assignBofO[i];
          int eid = assignAofO[e];
          // iresid = Btop.resids[iid];
          // eresid = Atop.resids[eid];
          iname = Btop.names[iid];
          ename = Atop.names[eid];
        }else if(assignAofO[e] == -1 &&
                 assignBofO[i] == -1) {
          int iid = assignAofO[i];
          int eid = assignBofO[e];
          // iresid = Atop.resids[iid];
          // eresid = Btop.resids[eid];
          iname = Atop.names[iid];
          ename = Btop.names[eid];
        }
        Ofs << "\t" << (e + 1);
      }
      Ofs << endl;
    }
    Ofs << endl;
  }

  if(!restraints_atomname.empty()) {
    Ofs << "#ifdef RESTR" << endl;
    Ofs << "[ position_restraints ]" << endl;
    set<string> restcond;
    size_t pos = 0;
    size_t pend;
    while(true) {
      string token;
      pend = restraints_atomname.find(',', pos);
      // here pend might be either position or npos
      if(pend == std::string::npos) {
        token = std::string(restraints_atomname, pos);
      }else{
        token = std::string(restraints_atomname, pos, pend - pos);
      }
      restcond.insert(token);
      if(pend == std::string::npos) {
        break;
      }else{
        pos = pend + 1;
      }
    }
    for(int i = 0; i < N; i++) {
      if((assignAofO[i] != -1 && restcond.count(Atop.names[assignAofO[i]]) > 0) ||
         (assignBofO[i] != -1 && restcond.count(Btop.names[assignBofO[i]]) > 0)) {
        int functype = 1;
        double forceconst = 1000.;
        Ofs << (i + 1) << " " << functype << " ";
        Ofs << forceconst << " " << forceconst << " " << forceconst << endl;
      }
    }
    Ofs << "#endif" << endl;
  }

  // Generate footer, way to go!
  Ofs << "[ system ]" << endl;
  Ofs << "Merged structure from " 
      << title
      << endl;
  Ofs << endl;

  Ofs << "[ molecules ]" << endl;
  Ofs << "merged 1" << endl;
  Ofs << endl;

}

static string find_file_at_exe_path(const string& exepath, const string& fname)
{
  size_t p = exepath.find_last_of('/');
  if(p == string::npos) {
    throw runtime_error(string("Failed to find ") + fname);
  }
  string exedir_with_slash = exepath.substr(0, p + 1);
  return exedir_with_slash + fname;
}

static vector<assigner_conditions> load_dictionary(const string& exepath, const string& fname)
{
  ifstream ifs(fname);
  if(!ifs) {
    ifs.close();
    ifs.clear();
    ifs.open(find_file_at_exe_path(exepath, fname));
    if(!ifs) {
      cerr << "Failed to find \"" << fname << "\"" << endl;
      throw runtime_error("load_dictionary");
    }
  }
  string line;
  vector<assigner_conditions> ret;
  while(getline(ifs, line)) {
    if(line.empty() || (line.length() >= 1 && line[0] == '#')) {
      continue;
    }
    if(line.find_first_not_of(" \t\n") == string::npos) {
      // blank
      continue;
    }
    istringstream is(line);
    string actionstr, resname1, atomname1, resname2, atomname2;
    is >> actionstr >> resname1 >> atomname1 >> resname2 >> atomname2;
    if(is.fail()) {
      cerr << "Failed to parse: \"" << line << "\"" << endl;
      throw runtime_error("load_forbid_assign");
    }
    assigner_action action = ASSIGN_NR;
    if(actionstr == "ACCEPT") {
      action = ASSIGN_ACCEPT;
    }else if(actionstr == "REJECT") {
      action = ASSIGN_REJECT;
    }
    ret.emplace_back(assigner_conditions{action, resname1, atomname1, resname2, atomname2});
  }
  return ret;
}

int main(int argc, char* argv[])
{
  cmdline::parser p;
  
  p.add("help", 0, "Print this message");
  p.add("verbose", 'v', "Be more verbose");
  p.add("quiet", 'q', "Suppress unnecessary information");
  p.add<double>("maxdist", 0, "Maximum distances to fit", false, 1.0);
  p.add<double>("force-bond", 0, "Cluster bond force constants", false, 5e+5); // Typical bond force constants: 2-5e5. 5e5 is for double bonds.
  p.add<double>("force-angular", 0, "Cluster rotation fixing force constants", false, 5e+2); // Typical angle fc: 3-7e2. 
  p.add<double>("force-dihedral", 0, "Cluster rotation fixing force constants", false, 1e+2); // Typical dihedral fc: varies, up to 5e1 for peptide bond dihed rotation
  p.add<int>("max-warning", 0, "Maximum number of allowed errors", false, 0);
  p.add("debug", 0, "Debug");

  p.add<string>("structureA", 'A', "PDB structure A", true);
  p.add<string>("structureB", 'B', "PDB structure B", true);
  p.add<string>("topologyA", 'a', ".top A", true);
  p.add<string>("topologyB", 'b', ".top B", true);
  p.add<string>("structureO", 'O', "PDB structure to output", true);
  p.add<string>("structureOA", 0, "PDB structure to output (optional, for stateA)", false);
  p.add<string>("structureOB", 0, "PDB structure to output (optional, for stateB)", false);
  p.add<string>("topologyO", 'o', ".top output", true);
  p.add<string>("assign-dictionary", 0, "Specify assignment matching pair", false, "assign-dictionary.txt");
  p.add<string>("generate-restraint", 0, "Comma-separated list. If non-empty, restraints are generated for the specified atom name, with 1000 kJ/mol/nm^2", false, "");
  p.add("protein", 0, "Use protein atom-matching instead of nucleic acids");
  p.add("no-best-fitting", 0, "Stop best-fitting two structures at the beginning of the program");
  p.add("connectivity", 0, "match atoms by connectivity");
  p.add("assign-by-name", 0, "match atoms by both connectivity and name");
  p.add("gen-exclusion", 0, "Program writes exclusions explicitly instead of nexcl");
  p.add("honor-resnames", 0, "Do not check structure if residue id & residue name matches, and use dictionary matching procedure");
  p.add("wang-approximation", 0, "Enable an approximation described in Supp. Info in Wang et al. 10.1073/pnas.1114017109");
  p.add("disable-cmap-error", 0, "FEP of CMAP requires a specially patched GROMACS. Only if you know what you are doing specify this flag.");
  p.add("disable-pairs-error", 0, "FEP of pairs section requires a specially patched GROMACS. Only if you know what you are doing specify this flag.");
  p.add("use-topological-min", 0, "To calculate the restraint bond length, angle and dihedral, use information in topology file rather than the structure. "
    "Useful when one of input structures are suboptimal.");

  bool ok = p.parse(argc, argv);

  if (!ok || p.exist("help")) {
    cerr << p.usage();
    return ((ok && p.exist("help")) ? 0 : 1);
  }
  verbose = p.exist("verbose");
  quiet = p.exist("quiet");

  pdb Apdb(p.get<string>("structureA"));
  pdb Bpdb(p.get<string>("structureB"));
  topology Atop(p.get<string>("topologyA"));
  topology Btop(p.get<string>("topologyB"));

  assert(Atop.names.size() == Apdb.get_atomnames().size());
  assert(Btop.names.size() == Bpdb.get_atomnames().size());

  Matrix3Xd Acoords, Bcoords;
  Acoords = Apdb.get_coords();
  Bcoords = Bpdb.get_coords();

  const vector<string>& Anames = Apdb.get_atomnames();
  const vector<string>& Bnames = Bpdb.get_atomnames();

  // fit Bpdb structure into Apdb by comparing CA or P-coordinates
  if(p.exist("protein")){
    VectorXd Amass = set_selected_mass(Anames, "CA");
    VectorXd Bmass = set_selected_mass(Bnames, "CA");

    fit_selected(Amass, Bmass, Acoords, Bcoords);
  }else if(p.exist("no-best-fitting")){
    // do nothing
  }else{
    VectorXd Amass = set_selected_mass(Anames, "P");
    VectorXd Bmass = set_selected_mass(Bnames, "P");

    fit_selected(Amass, Bmass, Acoords, Bcoords);
  }

  // make a distance matrix
  MatrixXd distmat(Acoords.cols(), Bcoords.cols());

  for(int ia = 0; ia < Acoords.cols(); ++ia) {
    char chainA = Apdb.get_chains()[ia];
    for(int jb = 0; jb < Bcoords.cols(); ++jb) {
      char chainB = Bpdb.get_chains()[jb];
      // Quick hack: prevent different chain atoms to be assigned the same atom
      distmat(ia, jb) = 
        (chainA == chainB ? 
         (Acoords.col(ia) - Bcoords.col(jb)).norm()
         : numeric_limits<double>::max());
    }
  }

  // assign atoms by distance
  vector<int> assignBofA(Acoords.cols(), -1);
  vector<int> assignAofB(Bcoords.cols(), -1);
  vector<int> Adepth;

  assigner_dictionary dictionary = load_dictionary(string(argv[0]), p.get<string>("assign-dictionary"));
  if(verbose) {
    cerr << "DEBUG: dictionary size: " << dictionary.size() << endl;
  }

  bool use_connectivity = p.exist("connectivity") || p.exist("assign-by-name");
  double pdist = p.get<double>("maxdist");
  int assigned;
  if(p.exist("honor-resnames")) {
    assign_atoms_resinfo(Atop, Btop, dictionary,
                         &assignBofA, &assignAofB);
  }else{
    if(p.exist("protein")) {
      assigned = assign_atoms("",    "CA",    Atop, Btop,
                              distmat, assignBofA, assignAofB, pdist);
    }else{
      assigned = assign_atoms("P",   nullptr, Atop, Btop,
                              distmat, assignBofA, assignAofB, pdist);
    }
    if(verbose) {
      cout << "Assigned " << assigned << " atoms @ mainchain phase" << endl;
    }
    if(assigned == 0) {
      cout << "Too few assigned atoms at mainchain phase" << endl;
      exit(1);
    }

    if(use_connectivity) {
      assign_atoms_connectivity(Atop, Btop, distmat, assignBofA, assignAofB, Adepth, pdist, p.exist("assign-by-name"));
      correct_assign_by_exclusion(Atop, Btop, Apdb, Bpdb, assignBofA, assignAofB, Adepth);
    }else{
      assign_atoms("NCO",  nullptr, Atop, Btop, distmat, assignBofA, assignAofB, pdist);
      assign_atoms("H",    nullptr, Atop, Btop, distmat, assignBofA, assignAofB, pdist);
      assign_atoms("HNCO", nullptr, Atop, Btop, distmat, assignBofA, assignAofB, pdist); // re-asign unassigned
      throw runtime_error("FIXME: here we need to check exclusion and exit on error");
    }
  }
  if(!quiet){
    cout << "Atoms assigned:" << endl;
    cout << "Assigned:" << endl;
    for(int i = 0; i < (int)Anames.size(); ++i) {
      if(assignBofA[i] != -1) {
        int j = assignBofA[i];
        if(!verbose &&
           Apdb.get_residuenames()[i] == Bpdb.get_residuenames()[j] &&
           Apdb.get_resids()[i] == Bpdb.get_resids()[j] &&
           Anames[i] == Bnames[j])
          continue; // suppress printing
        cout << " " 
             << Apdb.get_chains()[i] << ":"
             << Apdb.get_residuenames()[i] << ":"
             << Apdb.get_resids()[i] << ":" << Anames[i] 
             << " -> " 
             << Bpdb.get_chains()[j] << ":"
             << Bpdb.get_residuenames()[j] << ":"
             << Bpdb.get_resids()[j] << ":" << Bnames[j]
             << " dist = " << (Acoords.col(i) - Bcoords.col(j)).norm()
             << endl;
      }
    }
    cout << "Unssigned in A:" << endl;
    for(int i = 0; i < (int)Anames.size(); ++i) {
      if(assignBofA[i] == -1) {
        cout << " "
             << Apdb.get_chains()[i] << ":"
             << Apdb.get_residuenames()[i] << ":"
             << Apdb.get_resids()[i] << ":" << Anames[i] 
             << endl;
      }
    }
    cout << "Unassigned in B:" << endl;
    for(int i = 0; i < (int)Bnames.size(); ++i) {
      if(assignAofB[i] == -1) {
        cout << " " 
             << Bpdb.get_chains()[i] << ":"
             << Bpdb.get_residuenames()[i] << ":"
             << Bpdb.get_resids()[i] << ":" << Bnames[i] 
             << endl;
      }
    }
  }

  vector<int> assignOofA, assignOofB, assignAofO, assignBofO;
  int N;

  {
    // Make map from A/B, to/from output
    // First make residue index mapping
    vector<int> Aresindex, Bresindex;
    int NAres, NBres;
    make_residue_index(Apdb, &Aresindex, &NAres);
    make_residue_index(Bpdb, &Bresindex, &NBres);
    merge_assigns(assignBofA, assignAofB, 
                  Aresindex, Bresindex,
                  assignOofA, assignOofB,
                  assignAofO, assignBofO,
                  &N);
  }

  // Sanity check
  {
    map<int, set<int> > Aconnectivity, Bconnectivity;
    for(const auto &Abond: Atop.bonds) {
      int Aa = assignOofA[std::get<0>(Abond.first)];
      int Ab = assignOofA[std::get<1>(Abond.first)];
      Aconnectivity[Aa].insert(Ab);
      Aconnectivity[Ab].insert(Aa);
    }
    
    for(const auto &Bbond: Btop.bonds) {
      int Ba = assignOofB[std::get<0>(Bbond.first)];
      int Bb = assignOofB[std::get<1>(Bbond.first)];
      Bconnectivity[Ba].insert(Bb);
      Bconnectivity[Bb].insert(Ba);
    }

    int errcnt = 0;
    for(int i = 0; i < N; ++i) {
      const set<int> &As = Aconnectivity[i];
      const set<int> &Bs = Bconnectivity[i];
      if(!(std::includes(As.begin(), As.end(),
                         Bs.begin(), Bs.end()) ||
           std::includes(Bs.begin(), Bs.end(),
                         As.begin(), As.end())) &&
         !use_connectivity) {
        // FIXME TODO: this part really should use A/Bconn without dummy atoms.
        bool unmatched_are_all_dummy = true;
        for(int ao: As) {
          int b = assignBofO[ao];
          if(b == -1) {
            // A real B dummy
            continue;
          }
          if(Bs.count(assignOofB[b]) == 0) {
            // A real and B real, but does not appear in B connectivity
            unmatched_are_all_dummy = false;
          }
        }
        for(int bo: Bs) {
          int a = assignAofO[bo];
          if(a == -1) {
            // A dummy B real
            continue;
          }
          if(As.count(assignOofA[a]) == 0) {
            // A real and B real, but does not appear in A connectivity
            unmatched_are_all_dummy = false;
          }
        }
        if(unmatched_are_all_dummy) {
          // OK, ignore this
          continue;
        }
        errcnt++;
        int Ai = assignAofO[i];
        int Bi = assignBofO[i];
        if(Ai == -1 || Bi == -1) {
          cerr << "Unknown error: Ai / Bi = -1 in sanity check" << endl;
          cerr << "Ai = " << Ai << ", Bi = " << Bi << ", i = " << i << endl;
          continue;
        }
        cerr << "*** Warning: unmatched connectivity ***" << endl;
        cerr << Atop.resids[Ai] << ":" << Atop.names[Ai] << "(A)" << ": ";
        for(const auto &e: As) {
          int Ae = assignAofO[e];
          cerr << (Ae == -1 ? "PHA" : Atop.names[Ae]) << " ";
        }
        cerr << endl;
        cerr << Btop.resids[Bi] << ":" << Btop.names[Bi] << "(B)" << ": ";
        for(const auto &e: Bs) {
          int Be = assignBofO[e];
          cerr << (Be == -1 ? "PHA" : Btop.names[Be]) << " ";
        }
        cerr << endl;

        cerr << "Sanity check: " << endl;
        cerr << "i = " << i << ", Ai = " << Ai << ", Bi = " << Bi << endl;
        cerr << "BofA[Ai] = " << assignBofA[Ai] 
             << ", Anames[Ai] = " << Anames[Ai] << endl;
        cerr << "AofB[Bi] = " << assignAofB[Bi]
             << ", Bnames[Bi] = " << Bnames[Bi] << endl;
        cerr << "Raw As:";
        for(const auto &e: As) {
          cerr << " " << e;
        }
        cerr << endl;
        cerr << "Raw Bs:";
        for(const auto &e: Bs) {
          cerr << " " << e;
        }
        cerr << endl;
      }
    }
    if(errcnt > p.get<int>("max-warning")) {
      cerr << "Number of warning exceeds max-warning" << endl;
      exit(1);
    }

    if(Atop.cmaps.size() > 0 || Btop.cmaps.size() > 0) {
      if(!p.exist("disable-cmap-error")) {
        cerr << "FEP with CMAP requires modified GROMACS. If you have the patch and you understand that the modified version must be used, use --disable-cmap-error option." << endl;
        exit(1);
      }
    }
  }

  vector<double> stateA_weights(HAM_NR, 0.);
  vector<double> stateB_weights(HAM_NR, 1.);
  struct dummy_anchoring_parameters cparams;
  cparams.force_bond = p.get<double>("force-bond");
  cparams.force_angular = p.get<double>("force-angular");
  cparams.force_dihedral = p.get<double>("force-dihedral");

  vector<string> rest_excluded_atoms;
  string outfname = p.get<string>("topologyO");
  string title = p.get<string>("topologyA") + string(" and ") + p.get<string>("topologyB");
  vector<string> args(argv, argv + argc);
    
  generate_topology(outfname,
                    args,
                    title,
                    Atop,
                    Btop,
                    Acoords,
                    Bcoords,
                    stateA_weights,
                    stateB_weights,
                    assignOofA,
                    assignOofB, 
                    assignAofO,
                    assignBofO,
                    cparams,
                    1.0,
                    rest_excluded_atoms,
                    p.exist("gen-exclusion"),
                    p.exist("wang-approximation"),
                    p.exist("disable-pairs-error"),
                    p.exist("use-topological-min"),
                    p.get<string>("generate-restraint"));

  // Generate PDB
  enum { OPDB_MERGED, OPDB_A, OPDB_B, OPDB_END };
  const char *options[OPDB_END] = { "structureO", "structureOA", "structureOB" };
  for(int opdbtype = 0; opdbtype < OPDB_END; opdbtype++) {
    if(!p.exist(options[opdbtype])) continue;
    ofstream crdfs(p.get<string>(options[opdbtype]));
    crdfs.setf(ios::fixed);

    for(int i = 0; i < N; ++i) {
      Vector3d crd;
      int resid = 0;
      string atom, resname;

      int acrd = assignAofO[i];
      int bcrd = assignBofO[i];
      switch(opdbtype) {
      case OPDB_MERGED:
        if(assignAofO[i] != -1) {
          crd = Acoords.col(acrd);
          atom = Atop.names[acrd];
          resid = Atop.resids[acrd];
          resname = Atop.resnames[acrd];
        }else{
          assert(bcrd >= 0);
          crd = Bcoords.col(bcrd);
          atom = Btop.names[bcrd];
          resid = Btop.resids[bcrd];
          resname = Btop.resnames[bcrd];
        }
        break;
      case OPDB_A:
        if(assignAofO[i] != -1) {
          crd = Acoords.col(acrd);
          atom = Atop.names[acrd];
          resid = Atop.resids[acrd];
          resname = Atop.resnames[acrd];
        }else{
          assert(bcrd >= 0);
          crd = Bcoords.col(bcrd);
          atom = "DU";
          resid = Btop.resids[bcrd];
          resname = "DUM";
        }
        break;
      case OPDB_B:
        if(assignBofO[i] != -1) {
          crd = Bcoords.col(bcrd);
          atom = Btop.names[bcrd];
          resid = Btop.resids[bcrd];
          resname = Btop.resnames[bcrd];
        }else{
          assert(acrd >= 0);
          crd = Acoords.col(acrd);
          atom = "DU";
          resid = Atop.resids[acrd];
          resname = "DUM";
        }
        break;
      default:
        // do nothing
        break;
      }

      crdfs << setw(6) << "ATOM  " 
            << setw(5) << std::right << (i + 1)
            << " ";
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

  return 0;
}

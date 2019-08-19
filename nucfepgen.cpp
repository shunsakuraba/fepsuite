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

void output_dummies(ofstream &Ofs, const Matrix3Xd &Acoords, 
                    double force_bond,
                    double force_angluar,
                    double force_dihedral,
                    const vector<int>& assignAofO,
                    const vector<int>& assignBofO,
                    const vector<int>& assignOofA,
                    const topology &Atop, 
                    bool isA)
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

  Ofs << "[ bonds ]" << endl;
  
  for(int p: Aphantoms) {
    int pa = assignAofO[p];
    int bond_o = bond_atom[p];
    int ba = assignAofO[bond_o];

    double d = (Acoords.col(pa) - Acoords.col(ba)).norm();
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

  Ofs << "[ angles ]" << endl;

  for(int p: Aphantoms) {
    int pa = assignAofO[p];
    int bond_o = bond_atom[p];
    int ba = assignAofO[bond_o];
    int angle_o = angle_atom[p];
    int aa = assignAofO[angle_o];

    double ang_rad = angle(Acoords.col(pa), Acoords.col(ba), Acoords.col(aa));

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

  Ofs << "[ dihedrals ]" << endl;

  for(int p: Aphantoms) {
    int pa = assignAofO[p];
    int bond_o = bond_atom[p];
    int ba = assignAofO[bond_o];
    int angle_o = angle_atom[p];
    int aa = assignAofO[angle_o];
    int dihed_o = dihedral_atom[p];
    int da = assignAofO[dihed_o];

    double dih_rad = dihedral(Acoords.col(pa), Acoords.col(ba), Acoords.col(aa), Acoords.col(da));
    
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

static void merge_assigns(const vector<int>& assignBofA,
                          const vector<int>& assignAofB,
                          vector<int>& assignOofA,
                          vector<int>& assignOofB,
                          vector<int>& assignAofO,
                          vector<int>& assignBofO,
                          int &N)
{
  int Bunassigned = 0;
  for(int i = 0; i < (int)assignAofB.size(); ++i) {
    if(assignAofB[i] == -1) Bunassigned++;
  }
  N = assignBofA.size() + Bunassigned;
  assignOofA = vector<int>(assignBofA.size(), -1);
  assignOofB = vector<int>(assignAofB.size(), -1);
  assignAofO = vector<int>(N, -1);
  assignBofO = vector<int>(N, -1);

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
                                 vector<int> &assignBofA,
                                 vector<int> &assignAofB,
                                 vector<int> &Adepth)
{
  vector<vector<int> > Apairs_A(assignBofA.size()), Bpairs_B(assignAofB.size());
  Atop.convert_pairs_to_adj_list(Apairs_A);
  Btop.convert_pairs_to_adj_list(Bpairs_B);
  bool found;
  do {
    found = false;
    int maxdepth = -1;
    int Amax, Bmax = -1;
    vector<int> assignOofA, assignOofB, assignAofO, assignBofO;
    int N;
    // try to generate current assign
    merge_assigns(assignBofA, assignAofB, 
                  assignOofA, assignOofB,
                  assignAofO, assignBofO,
                  N);
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
           << Btop.names[Bmax] << endl;
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
                       bool gen_exclusion)
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
       Atop.defaults == topology::topotype::CHARMM)) {
    cerr << "Unsupported defaults type" << endl;
    exit(1);
  }
  Ofs << "[ defaults ]" << endl;
  Ofs << "; nbfunc        comb-rule       gen-pairs       fudgeLJ fudgeQQ" << endl;
  if(Atop.defaults == topology::topotype::AMBER) {
    Ofs << "1               2               yes             0.5     0.8333" << endl;
  }else if(Atop.defaults == topology::topotype::CHARMM) {
    Ofs << "1 2 yes 1.0 1.0" << endl;
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
    if(func != 1) {
      throw runtime_error("Bond func != 1 not supported");
    }
    const vector<double>& Afactors = v.second;
    vector<double> Bfactors(2, 0);
    if(Btop.bonds.count(keyB) > 0) {
      Bfactors = Btop.bonds.at(keyB);
    }else if(Afactors.size() > 0){
      Bfactors[0] = Afactors[0];
      Bfactors[1] = 0.0;
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
          << " " << 0.00;
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
            << " " << setw(12) << 0.00;
      }else{
        Ofs << " " << setw(12) << Bfactors[0]
            << " " << setw(12) << 0.00
            << " " << setw(12) << Bfactors[2]
            << " " << setw(12) << 0.00;
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

  // pairs
  Ofs << "[ pairs ]" << endl;
  for(const auto& v: Atop.pairs) {
    auto k = v.first;
    int x = k.first;
    int y = k.second;
    int func = v.second;
    
    int x_in_O = assignOofA[x];
    int y_in_O = assignOofA[y];
    if(func != 1) {
      throw runtime_error("pairs func != 1 not supported");
    }
    Ofs << x_in_O + 1 << " " << y_in_O + 1 << " " 
        << func << endl;
  }
  // Remaining: not exists in A but exists in B
  for(const auto& v: Btop.pairs) {
    auto k = v.first;
    int x = k.first;
    int y = k.second;
    int func = v.second;
    
    int x_in_O = assignOofB[x];
    int y_in_O = assignOofB[y];
    auto keyA = 
      make_pair(assignAofO[x_in_O],
                assignAofO[y_in_O]);
    if(Atop.pairs.count(keyA) > 0) {
      continue;
    }
    Ofs << x_in_O + 1 << " " << y_in_O + 1 << " " 
        << func << endl;
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
                   Atop, true);
    Ofs << "; Atoms existing in B state" << endl;
    output_dummies(Ofs, Bcoords, 
                   cparams.force_bond,
                   cparams.force_angular,
                   cparams.force_dihedral,
                   assignBofO,
                   assignAofO,
                   assignOofB,
                   Btop, false);
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
  p.add("protein", 0, "Use protein atom-matching instead of nucleic acids");
  p.add("connectivity", 0, "match atoms by connectivity");
  p.add("assign-by-name", 0, "match atoms by both connectivity and name");
  p.add("gen-exclusion", 0, "Program writes exclusions explicitly instead of nexcl");
  p.add("honor-resnames", 0, "Do not check structure if residue id & residue name matches");

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

  // fit Bpdb structure into Apdb by comparing P-coordinates
  if(p.exist("protein")){
    VectorXd Amass = set_selected_mass(Anames, "CA");
    VectorXd Bmass = set_selected_mass(Bnames, "CA");

    fit_selected(Amass, Bmass, Acoords, Bcoords);
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

  bool use_connectivity = p.exist("connectivity") || p.exist("assign-by-name");
  double pdist = p.get<double>("maxdist");
  int assigned;
  if(p.exist("protein")) {
    assigned = assign_atoms("",    "CA",    Anames, Bnames, distmat, assignBofA, assignAofB, pdist);
  }else{
    assigned = assign_atoms("P",   nullptr, Anames, Bnames, distmat, assignBofA, assignAofB, pdist);
  }
  if(verbose) {
    cout << "Assigned " << assigned << " atoms @ mainchain phase" << endl;
  }
  if(assigned == 0) {
    cout << "Too few assigned atoms at mainchain phase" << endl;
    exit(1);
  }
  if(p.exist("honor-resnames")) {
    assign_atoms_resinfo(Anames, Bnames, Apdb.get_residuenames(), Bpdb.get_residuenames(),
                         Apdb.get_resids(), Bpdb.get_resids(),
                         assignBofA, assignAofB);
  }

  if(use_connectivity) {
    assign_atoms_connectivity(distmat, Atop, Btop, assignBofA, assignAofB, Adepth, pdist, p.exist("assign-by-name"));
    correct_assign_by_exclusion(Atop, Btop, assignBofA, assignAofB, Adepth);
  }else{
    assign_atoms("NCO",  nullptr, Anames, Bnames, distmat, assignBofA, assignAofB, pdist);
    assign_atoms("H",    nullptr, Anames, Bnames, distmat, assignBofA, assignAofB, pdist);
    assign_atoms("HNCO", nullptr, Anames, Bnames, distmat, assignBofA, assignAofB, pdist); // re-asign unassigned
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
           Anames[i] == Anames[j])
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

  // Make map from A/B, to/from output
  merge_assigns(assignBofA, assignAofB, 
                assignOofA, assignOofB,
                assignAofO, assignBofO,
                N);


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
                    p.exist("gen-exclusion"));

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

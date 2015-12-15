#include "cmdline.h"
#include "bestfit.hpp"
#include "pdb.hpp"
#include "topology.hpp"
#include "assign.hpp"
#include "select.hpp"

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

bool verbose = false;
bool quiet = false;

void output_clusters(ofstream &Ofs, const Matrix3Xd &Acoords, 
                     double maxdist,
                     double force_connect,
                     double force_within_cluster,
                     const vector<int>& assignAofO,
                     const vector<int>& assignBofO,
                     const topology &Atop, 
                     bool isA)
{
  int N = (int)assignAofO.size();
  vector<int> Aphantoms; // A's atoms which is phantom in state B, numbered in O

  for(int i = 0; i < N; ++i) {
    if(assignBofO[i] == -1) {
      Aphantoms.push_back(i);
    }
  }

  // make clusters of phantom atoms
  disjoint_set As(N);
  for(int x: Aphantoms) {
    assert(assignAofO[x] >= 0);
    Vector3d apt = Acoords.col(assignAofO[x]);
    for(int y: Aphantoms) {
      assert(assignAofO[y] >= 0);
      Vector3d bpt = Acoords.col(assignAofO[y]);
      if(x == y) continue;
      if((apt - bpt).squaredNorm() > maxdist * maxdist) {
        continue;
      }
      if(Atop.resids[assignAofO[x]] != Atop.resids[assignAofO[y]]) {
        continue;
      }
      As.join(x, y);
    }
  }
  
  // Turn to clusters
  multimap<int, int> Aclusters;
  for(int x: Aphantoms) {
    const int p = As.find_parent(x);
    Aclusters.insert(make_pair(p, x));
  }

  // find the closest cluster-to-non-phantom-atom distance
  map<int, pair<int, double> > Aclusterdist;

  for(const auto &x: Aclusters) {
    Aclusterdist[x.first] = 
      make_pair(-1, numeric_limits<double>::max());
  }
    
  for(int x: Aphantoms) {
    const int p = As.find_parent(x);
    int curminfrom = Aclusterdist[p].first;
    double curmin  = Aclusterdist[p].second;

    Vector3d apt = Acoords.col(assignAofO[x]);
    for(int y = 0; y < N; ++y) {
      if(assignAofO[y] == -1 || assignBofO[y] == -1) {
        // Both A phantom and B phantom are excluded
        continue;
      }
      // different resid ==> rejected
      if(Atop.resids[assignAofO[y]] != Atop.resids[assignAofO[x]]) {
        continue;
      }
      Vector3d bpt = Acoords.col(assignAofO[y]);
      double d = (apt - bpt).norm();
      if(d < curmin) {
        curminfrom = y;
        curmin = d;
      }
    }
    Aclusterdist[p] = make_pair(curminfrom, curmin);
  }
  
  for(const auto& x: Aclusterdist) {
    const double scale = 0.1; // angstrom to nm
    int xa = assignAofO[x.second.first];
    int par = x.first;
    auto ar = Aclusters.equal_range(par);
    for(auto it = ar.first; it != ar.second; ++it) {
      {
        int ya = assignAofO[it->second];
        double d = (Acoords.col(xa) - Acoords.col(ya)).norm();
        Ofs << setw(6) << x.second.first + 1 << " "
            << setw(6) << it->second + 1 << " "
            << setw(4) << 6 << " ";
        if(isA) {
          Ofs << setw(8) << d * scale << " "
              << setw(8) << 0.0 << " ";
        }
        Ofs << setw(8) << d * scale << " "
            << setw(8) << force_connect << " ";
        if(!isA) {
          Ofs << setw(8) << d * scale << " "
              << setw(8) << 0.0 << " ";
        }        
        Ofs << "; (to-cluster) " 
            << setw(3) << Atop.resids[xa] << " "
            << setw(5) << Atop.names[xa] << " - "
            << setw(3) << Atop.resids[ya] << " "
            << setw(5) << Atop.names[ya]
            << endl;
      }
      
      for(auto it2 = ar.first; it2 != ar.second; ++it2) {
        if(it->second >= it2->second) continue;
        int xa = assignAofO[it->second];
        int ya = assignAofO[it2->second];
        double d = (Acoords.col(xa) - 
                    Acoords.col(ya)).norm();
        
        Ofs << setw(6) << it->second + 1 << " "
            << setw(6) << it2->second + 1 << " "
            << setw(4) << 6 << " ";
        if(isA) {
          Ofs << setw(8) << d * scale << " "
              << setw(8) << 0.0 << " ";
        }
        Ofs << setw(8) << d * scale << " "
            << setw(8) << force_within_cluster << " ";
        if(!isA) {
          Ofs << setw(8) << d * scale << " "
              << setw(8) << 0.0 << " ";
        }
        Ofs << "; (in-cluster) " 
            << setw(3) << Atop.resids[xa] << " "
            << setw(5) << Atop.names[xa] << " - "
            << setw(3) << Atop.resids[ya] << " "
            << setw(5) << Atop.names[ya]
            << endl;
      }
    }
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
    int Amax, Bmax;
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

int main(int argc, char* argv[])
{
  cmdline::parser p;
  
  p.add("help", 0, "Print this message");
  p.add("verbose", 'v', "Be more verbose");
  p.add("quiet", 'q', "Suppress unnecessary information");
  p.add<double>("maxdist", 0, "Maximum distances to fit", false, 1.0);
  p.add<double>("cluster-dist", 0, "Maximum distances for phantom atoms clustering", false, 2.0);
  p.add<double>("force-constant", 0, "Cluster fixing force constants", false, 5e+4);
  p.add<double>("force-within-cluster", 0, "Cluster force constants", false, 5e+4);
  p.add<int>("max-warning", 0, "Maximum number of allowed errors", false, 0);
  p.add("debug", 0, "Debug");

  p.add<string>("structureA", 'A', "PDB structure A", true);
  p.add<string>("structureB", 'B', "PDB structure B", true);
  p.add<string>("topologyA", 'a', ".top A", true);
  p.add<string>("topologyB", 'b', ".top B", true);
  p.add<string>("structureO", 'O', "PDB structure to output", true);
  p.add<string>("topologyO", 'o', ".top output", true);
  p.add("connectivity", 0, "match atoms by connectivity");
  p.add("assign-by-name", 0, "match atoms by both connectivity and name");
  p.add("gen-exclusion", 0, "Program writes exclusions explicitly instead of nexcl");

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

  Matrix3Xd Acoords, Bcoords;
  Acoords = Apdb.get_coords();
  Bcoords = Bpdb.get_coords();

  const vector<string>& Anames = Apdb.get_atomnames();
  const vector<string>& Bnames = Bpdb.get_atomnames();

  // fit Bpdb structure into Apdb by comparing P-coordinates
  {
    
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
  assign_atoms("P",   Anames, Bnames, distmat, assignBofA, assignAofB, pdist);
  if(use_connectivity) {
    assign_atoms_connectivity(distmat, Atop, Btop, assignBofA, assignAofB, Adepth, pdist, p.exist("assign-by-name"));
    correct_assign_by_exclusion(Atop, Btop, assignBofA, assignAofB, Adepth);
  }else{
    assign_atoms("NCO", Anames, Bnames, distmat, assignBofA, assignAofB, pdist);
    assign_atoms("H",   Anames, Bnames, distmat, assignBofA, assignAofB, pdist);
    assign_atoms("HNCO", Anames, Bnames, distmat, assignBofA, assignAofB, pdist); // re-asign unassigned
  }
  if(!quiet){
    cout << "Atoms assigned:" << endl;
    cout << "Assigned:" << endl;
    for(int i = 0; i < (int)Anames.size(); ++i) {
      if(assignBofA[i] != -1) {
        int j = assignBofA[i];
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
        cerr << Atop.names[Ai] << "(A)" << ": ";
        for(const auto &e: As) {
          int Ae = assignAofO[e];
          cerr << (Ae == -1 ? "PHA" : Atop.names[Ae]) << " ";
        }
        cerr << endl;
        cerr << Btop.names[Bi] << "(B)" << ": ";
        for(const auto &e: Bs) {
          int Be = assignBofO[e];
          cerr << (Be == -1 ? "PHA" : Btop.names[Be]) << " ";
        }
        cerr << endl;
      }
    }
    if(errcnt > p.get<int>("max-warning")) {
      cerr << "Number of warning exceeds max-warning" << endl;
      exit(1);
    }
  }

  // output 
  ofstream Ofs(p.get<string>("topologyO"));
  Ofs.setf(ios::scientific);
  Ofs.setf(ios::right);
  Ofs << setprecision(5);

  Ofs << "; generated by:" << endl;
  Ofs << ";  " << argv[0];
  for(int i = 1; i < argc; ++i) {
    Ofs << " " << argv[i];
  }
  Ofs << endl;
  
  // Defaults section
  if(Atop.defaults != topology::topotype::AMBER ||
     Btop.defaults != topology::topotype::AMBER) {
    cerr << "Unsupported defaults type" << endl;
    exit(1);
  }
  Ofs << "[ defaults ]" << endl;
  Ofs << "; nbfunc        comb-rule       gen-pairs       fudgeLJ fudgeQQ" << endl;
  Ofs << "1               2               yes             0.5     0.8333" << endl;
  Ofs << endl;

  // atomtypes section
  Ofs << "[ atomtypes ]" << endl;
  for(const auto& v: Atop.atomtypes) {
    Ofs << v.second << endl;
  }
  for(const auto& v: Btop.atomtypes) {
    if(Atop.atomtypes.count(v.first) == 0) {
      Ofs << v.second << endl;
    }
  }
  // phantom type
  Ofs << "PHA PHA 0.0000 0.0000 A 0.0000 0.0000" << endl;
  Ofs << endl;

  // moleculetype section
  Ofs << "[ moleculetype ]" << endl;
  Ofs << ";name  nrexcl" << endl;
  {
    int nexcl = Atop.nexcl;
    if(p.exist("gen-exclusion")) {
      nexcl = 0;
    }
    Ofs << "merged " << nexcl << endl;
  }
  Ofs << endl;

  // atoms section
  Ofs << "[ atoms ]" << endl;
  Ofs << ";   nr  type  resi  res  atom  cgnr     charge      mass       typeB chargeB massB" << endl;

  for(int i = 0; i < N; ++i) {
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
          << " " << setw(12) << Atop.masses[aptr];
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
          << " " << setw(12) << 12.000;
    }

    // Fill in B-state
    if(assignBofO[i] != -1) {
      int bptr = assignBofO[i];
      Ofs << " " << setw(3) << Btop.types[bptr]
          << " " << setw(12) << Btop.charges[bptr]
          << " " << setw(12) << Btop.masses[bptr];
    }else{
      Ofs << " " << "PHA"
          << " " << setw(12) << 0.000
          << " " << setw(12) << 12.000;
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
      Bfactors = Btop.bonds[keyB];
    }else{
      Bfactors[0] = Afactors[0];
      Bfactors[1] = 0.0;
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
        << func
        << " " << Bfactors[0]
        << " " << 0.00;
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
    if(func != 1) {
      throw runtime_error("Angle func != 1 not supported");
    }
    const vector<double>& Afactors = v.second;
    vector<double> Bfactors(2, 0);
    if(Btop.angles.count(keyB) > 0) {
      Bfactors = Btop.angles[keyB];
    }else{
      Bfactors[0] = Afactors[0];
      Bfactors[1] = 0.0;
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
        << func
        << " " << setw(12) << Bfactors[0]
        << " " << setw(12) << 0.00;
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
    if(func != 1 && func != 3) {
      throw runtime_error("dihed func not in {1, 3} not supported");
    }
    const vector<double>& Afactors = v.second;
    vector<double> Bfactors(Afactors);
    if(Btop.diheds.count(keyB) > 0) {
      Bfactors = Btop.diheds[keyB];
    }else{
      for(auto &&v: Bfactors) {
        v = 0.0;
      }
      if(func == 1) {
        Bfactors[0] = Afactors[0];
        Bfactors[2] = Afactors[2];
      }
    }
    // output A-listed angles first
    Ofs << x_in_O + 1 << " " << y_in_O + 1 << " " 
        << z_in_O + 1 << " " << w_in_O + 1 << " "
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
    if(func == 1) {
      Ofs << " " << setw(12) << Bfactors[0]
          << " " << setw(12) << 0.00
          << " " << setw(12) << Bfactors[2];
    }else{
      for(double v: Bfactors) {
        (void) v;
        Ofs << " " << setw(12) << 0.00;
      }
    }
    for(auto v: Bfactors) {
      Ofs << " " << setw(12) << v;
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

  // bonds for phantom atoms
  Ofs << "[ bonds ]" << endl;
  {
    Ofs << "; Atoms existing in A state" << endl;
    output_clusters(Ofs, Acoords, 
                    p.get<double>("cluster-dist"),
                    p.get<double>("force-constant"),
                    p.get<double>("force-within-cluster"),
                    assignAofO,
                    assignBofO,
                    Atop, true);
    Ofs << "; Atoms existing in B state" << endl;
    output_clusters(Ofs, Bcoords, 
                    p.get<double>("cluster-dist"),
                    p.get<double>("force-constant"),
                    p.get<double>("force-within-cluster"),
                    assignBofO,
                    assignAofO,
                    Btop, false);
  }
  Ofs << endl;

  if(p.exist("gen-exclusion")) {
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
        int iresid = -1, eresid = -1;
        string iname, ename;
        if(assignAofO[i] == -1 &&
           assignBofO[e] == -1) {
          int iid = assignBofO[i];
          int eid = assignAofO[e];
          iresid = Btop.resids[iid];
          eresid = Atop.resids[eid];
          iname = Btop.names[iid];
          ename = Atop.names[eid];
        }else if(assignAofO[e] == -1 &&
                 assignBofO[i] == -1) {
          int iid = assignAofO[i];
          int eid = assignBofO[e];
          iresid = Atop.resids[iid];
          eresid = Btop.resids[eid];
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
      << p.get<string>("topologyA") 
      << " and " 
      << p.get<string>("topologyB")
      << endl;
  Ofs << endl;

  Ofs << "[ molecules ]" << endl;
  Ofs << "merged 1" << endl;
  Ofs << endl;

  // Generate PDB
  {
    ofstream crdfs(p.get<string>("structureO"));
    crdfs.setf(ios::fixed);

    for(int i = 0; i < N; ++i) {
      Vector3d crd;
      int resid;
      string atom, resname;
      
      if(assignAofO[i] != -1) {
        int acrd = assignAofO[i];
        crd = Acoords.col(acrd);
        atom = Atop.names[acrd];
        resid = Atop.resids[acrd];
        resname = Atop.resnames[acrd];
      }else{
        int bcrd = assignBofO[i];
        assert(bcrd >= 0);
        crd = Bcoords.col(bcrd);
        atom = Btop.names[bcrd];
        resid = Btop.resids[bcrd];
        resname = Btop.resnames[bcrd];
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

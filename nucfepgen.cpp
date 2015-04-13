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
        if(!isA) {
          Ofs << setw(8) << d << " "
              << setw(8) << 0.0 << " ";
        }
        Ofs << setw(8) << d << " "
            << setw(8) << force_connect << " ";
        if(isA) {
          Ofs << setw(8) << d << " "
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
        if(!isA) {
          Ofs << setw(8) << d << " "
              << setw(8) << 0.0 << " ";
        }
        Ofs << setw(8) << d << " "
            << setw(8) << force_within_cluster << " ";
        if(isA) {
          Ofs << setw(8) << d << " "
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

  double pdist = p.get<double>("maxdist");
  assign_atoms("P",   Anames, Bnames, distmat, assignBofA, assignAofB, pdist);
  assign_atoms("NCO", Anames, Bnames, distmat, assignBofA, assignAofB, pdist);
  assign_atoms("H",   Anames, Bnames, distmat, assignBofA, assignAofB, pdist);
  assign_atoms("HNCO", Anames, Bnames, distmat, assignBofA, assignAofB, pdist); // re-asign unassigned

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

  int Bunassigned = 0;
  for(int i = 0; i < (int)Bnames.size(); ++i) {
    if(assignAofB[i] == -1) Bunassigned++;
  }
  int N = Acoords.cols() + Bunassigned;
  
  // Make map from A/B, to/from output
  vector<int> assignOofA(Acoords.cols(), -1);
  vector<int> assignOofB(Bcoords.cols(), -1);
  vector<int> assignAofO(N, -1);
  vector<int> assignBofO(N, -1);

  for(int i = 0; i < (int)assignOofA.size(); ++i) {
    assignOofA[i] = i;
    assignAofO[i] = i;
  }
  for(int j = 0, ptr = Acoords.cols();
      j < (int)assignOofB.size(); ++j) {
    if(assignAofB[j] == -1) {
      assignOofB[j] = ptr;
      ptr++;
    }else{
      assignOofB[j] = assignOofA[assignAofB[j]];
    }
    assignBofO[assignOofB[j]] = j;
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
                         As.begin(), As.end()))) {
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
  Ofs << "merged " << Atop.nexcl << endl;
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
    
    int x_in_O = assignOofA[x];
    int y_in_O = assignOofA[y];
    int z_in_O = assignOofA[z];
    int w_in_O = assignOofA[w];
    topology::dihedkeytype keyB = 
      make_tuple(assignBofO[x_in_O],
                 assignBofO[y_in_O],
                 assignBofO[z_in_O],
                 assignBofO[w_in_O],
                 func);
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
    
    int x_in_O = assignOofB[x];
    int y_in_O = assignOofB[y];
    int z_in_O = assignOofB[z];
    int w_in_O = assignOofB[w];
    topology::dihedkeytype keyA = 
      make_tuple(assignAofO[x_in_O],
                 assignAofO[y_in_O],
                 assignAofO[z_in_O],
                 assignAofO[w_in_O],
                 func);
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
    output_clusters(Ofs, Acoords, 
                    p.get<double>("cluster-dist"),
                    p.get<double>("force-constant"),
                    p.get<double>("force-within-cluster"),
                    assignAofO,
                    assignBofO,
                    Atop, true);
    output_clusters(Ofs, Bcoords, 
                    p.get<double>("cluster-dist"),
                    p.get<double>("force-constant"),
                    p.get<double>("force-within-cluster"),
                    assignBofO,
                    assignAofO,
                    Btop, false);
  }
  Ofs << endl;

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
            << " "
            << (atom.length() == 4 ? atom.substr(3, 1) : " ")
            << setw(3) << std::left
            << (atom.length() == 4 ? atom.substr(0, 3) : atom)
            << " " // altLoc
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

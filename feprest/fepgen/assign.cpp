#include "assign.hpp"
#include <functional>
#include <limits>
#include <queue>

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

bool dict_match(const assigner_dictionary& dictionary,
                const string &Aresname,
                const string &Aname,
                const string &Bresname,
                const string &Bname)
{
  for(const assigner_conditions &ad: dictionary){
    for(int i = 0; i < 2; ++i) {
      const string& resname1 = (i == 0 ? Aresname : Bresname);
      const string& atomname1 = (i == 0 ? Aname : Bname);
      const string& resname2 = (i == 0 ? Bresname : Aresname);
      const string& atomname2 = (i == 0 ? Bname : Aname);
      
      if((ad.res1 == "*" || ad.res1 == resname1) &&
        (ad.res2 == "*" || ad.res2 == resname2) &&
        (ad.atom1 == "*" || ad.atom1 == atomname1) &&
        (ad.atom2 == "*" || ad.atom2 == atomname2)) {
        assigner_action action = ad.action;
        if(action == ASSIGN_ACCEPT) {
          return true;
        }
        if(action == ASSIGN_REJECT) {
          return false;
        }
      }
    }
  }
  return false;
}

int assign_atoms(const string& process_atoms,
                 const char* atomname,
                 const topology &Atop,
                 const topology &Btop,
                 const MatrixXd& distmat,
                 vector<int>& assignBofA,
                 vector<int>& assignAofB,
                 double threshold)
{
  const vector<string> &Anames = Atop.names;
  const vector<string> &Bnames = Btop.names;

  vector<bool> Aenable(Anames.size());
  vector<bool> Benable(Bnames.size());

  int assigned = 0;
  for(int i = 0; i < (int)Anames.size(); i++) {
    char at = atomtype(Anames[i]);
    if(atomname) {
      Aenable[i] = (Anames[i] == atomname);
    }else{
      Aenable[i] = at && (process_atoms.find(at) != string::npos);
    }
  }

  for(int i = 0; i < (int)Bnames.size(); i++) {
    char at = atomtype(Bnames[i]);
    if(atomname) {
      Benable[i] = (Bnames[i] == atomname);
    }else{
      Benable[i] = at && (process_atoms.find(at) != string::npos);
    }
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
      ++assigned;
    }
  }
  return assigned;
}

void assign_atoms_connectivity(const topology &Atop,
                               const topology &Btop,
                               const MatrixXd& distmat,
                               vector<int>& assignBofA,
                               vector<int>& assignAofB,
                               vector<int>& depth,
                               double threshold,
                               bool must_be_identical_names)
{
  // Turn bonds into adjacent list
  vector<vector<int> > Aadjlist, Badjlist;
  Atop.convert_bonds_to_adj_list(Aadjlist);
  Btop.convert_bonds_to_adj_list(Badjlist);
  
  depth = vector<int>(Aadjlist.size(), -1);
  
  // um, indeed, this need not be priority queue (BFS suffice)
  priority_queue<pair<int, int>, 
                 vector<pair<int, int> >,
                 greater<pair<int, int> > > pq;
  for(int i = 0; i < (int) assignBofA.size(); i++) {
    if(assignBofA[i] != -1) {
      pq.emplace(make_pair(0, i));
    }
  }

  while(!pq.empty()) {
    pair<int, int> e = pq.top();
    pq.pop();
    int Ai = e.second;

    if(depth[Ai] >= 0) continue;
    depth[Ai] = e.first;

    assert(assignBofA[Ai] >= 0);
    int Bi = assignBofA[Ai];
    assert(assignAofB[Bi] >= 0);
    
    for(int An: Aadjlist[Ai]) {
      // An is a neighbor of Ai
      if(assignBofA[An] != -1) continue;

      // Find best maching Bn s.t.
      // Bn is a neighbor of Bi
      // Bn is the closest to An and closer than the threshold
      // Bn is unassigned
      double mindist = numeric_limits<double>::max();
      int Bn = -1;
      for(int Bn_cand: Badjlist[Bi]) {
        if(assignAofB[Bn_cand] != -1) {
          continue;
        }
        double d = distmat(An, Bn_cand);
        if(d < mindist && d < threshold) {
          Bn = Bn_cand;
          mindist = d;
        }
      }
      if(Bn != -1) {
        // found Bn
        assert(An < (int)Atop.names.size());
        assert(Bn < (int)Btop.names.size());
        assert(An >= 0);
        assert(Bn >= 0);
        if(must_be_identical_names && 
           Atop.names[An] != Btop.names[Bn]) {
          continue;
        }
        assignBofA[An] = Bn;
        assignAofB[Bn] = An;
        pq.emplace(make_pair(e.first + 1, An));
      }
    }
  }

}


int assign_atoms_resinfo(const topology &Atop,
                         const topology &Btop,
                         const assigner_dictionary& dictionary,
                         vector<int>* assignBofA,
                         vector<int>* assignAofB)
{
  const vector<string> &Anames = Atop.names;
  const vector<string> &Bnames = Btop.names;
  const vector<string> &Aresnames = Atop.resnames;
  const vector<string> &Bresnames = Btop.resnames;
  const vector<int> &Aresids = Atop.resids;
  const vector<int> &Bresids = Btop.resids;
  const vector<int> &Achainids = Atop.chainids;
  const vector<int> &Bchainids = Btop.chainids;
  int assigned = 0;

  // Residues are grouped based on chainid. This is necessary for multi-chain systems (see Issue #4)
  multimap<pair<int, int>, size_t> Bresids_to_Batoms;
  for(int j = 0; j < (int)Bnames.size(); ++j) {
    int chainid = Bchainids[j];
    int resid = Bresids[j];
    Bresids_to_Batoms.insert({make_pair(chainid, resid), (size_t)j});
  }

  for(int i = 0; i < (int)Anames.size(); ++i) {
    if((*assignBofA)[i] != -1) {
      continue;
    }
    int found = -1;
    auto range = Bresids_to_Batoms.equal_range(make_pair(Achainids[i], Aresids[i]));
    for(auto p = range.first; p != range.second; ++p) {
      int j = p->second;
      if((*assignAofB)[j] != -1) {
        continue;
      }

      if((Aresnames[i] == Bresnames[j]) &&
         (Anames[i] == Bnames[j])) {
        found = j;
        break;
      }
      if(Aresnames[i] != Bresnames[j]) {
        if(dict_match(dictionary, Aresnames[i], Anames[i], Bresnames[j], Bnames[j])) {
          found = j;
          break;
        }
      }
    }
    if(found >= 0) {
      (*assignBofA)[i] = found;
      (*assignAofB)[found] = i;
      ++assigned;
    }
  }
  return assigned;
}



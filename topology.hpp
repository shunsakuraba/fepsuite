#pragma once

#include <map>
#include <string>
#include <vector>
#include <tuple>

struct topology {
  enum topotype { AMBER, CHARMM, OPLS };

  // a, b, (c, d,) func (, addendum)
  // um, i know, it's ugly
  typedef std::tuple<int, int, int> bondkeytype;
  typedef std::tuple<int, int, int, int> anglekeytype;
  typedef std::tuple<int, int, int, int, int, int> dihedkeytype;
  typedef std::tuple<std::string, std::string, int> pairtype;
  typedef std::tuple<std::string, std::string, int> bondtype;
  typedef std::tuple<std::string, std::string, std::string, int> angletype;
  typedef std::tuple<std::string, std::string, std::string, std::string, int> dihedraltype;

  topotype defaults;
  std::string moleculename;
  int nexcl;

  // atomtypes
  // stores lines
  std::map<std::string, std::string> atomtypes;
  // pairtypes, stores the whole line
  std::map<pairtype, std::string> pairtypes;
  // bondtypes
  std::map<bondtype, std::vector<std::string>> bondtypes;
  // angletypes
  std::map<angletype, std::vector<std::string>> angletypes;
  // dihedraltypes
  std::map<dihedraltype, std::vector<std::string>> dihedraltypes;
  
  // atoms
  std::vector<std::string> names;
  std::vector<std::string> types;
  std::vector<int> resids;
  std::vector<std::string> resnames;
  std::vector<double> charges;
  std::vector<double> masses;

  // bonds. keytype atom numbers are 0-origin
  std::map<bondkeytype, std::vector<double>> bonds;
  
  // (1-4) pairs
  std::map<std::pair<int, int>, int> pairs;
  
  // angles
  std::map<anglekeytype, std::vector<double>> angles;

  // dihedrals
  std::map<dihedkeytype, std::vector<double>> diheds;


  topology();
  topology(const std::string& fname);
  
  void write(const std::string& fname);

  void convert_bonds_to_adj_list(std::vector<std::vector<int> > &adj_list) const;
  void convert_pairs_to_adj_list(std::vector<std::vector<int> > &pairs) const;
};


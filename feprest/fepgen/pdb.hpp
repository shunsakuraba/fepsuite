#pragma once

#include <Eigen/Core>
#include <vector>
#include <string>
#include <iosfwd>
#include <memory>

class pdb
{
public:
  pdb();

  /*
    Read pdb file. When altLoc is specified,
    it reads the line without altLoc and specified altLoc.
   */
  pdb(const std::string& fname, char altloc = '\0');

  const Eigen::Matrix3Xd& get_coords() const { return coords; }
  const Eigen::VectorXd& get_betas() const { return betas; }
  const Eigen::VectorXd& get_occupancies() const { return occupancies; }
  const std::vector<char>& get_chains() const { return chains; }
  const std::vector<std::string>& get_residuenames() const { return residuenames; }
  const std::vector<int>& get_resids() const { return resids; }
  const std::vector<std::string>& get_atomnames() const { return atomnames; }
  const size_t get_numatoms() const { return (size_t)coords.cols(); } 

  /**
     Try to read the next structure from the specified PDB.
     Returns false if no more structures are found.
     In this case, the structures in this class does not change.
     If a structure is found, returns true, and the system state is overwritten by the new coordinates.
     When PDB file is corrupt, its behaivour is undefined.
   */
  bool next();

private:

  std::shared_ptr<std::ifstream> fh;
  char altloc;
  
  Eigen::Matrix3Xd coords;
  Eigen::VectorXd betas;
  Eigen::VectorXd occupancies;

  std::vector<char> chains;
  std::vector<std::string> residuenames;
  std::vector<int> resids;
  std::vector<std::string> atomnames;

  bool load_pdb();
};



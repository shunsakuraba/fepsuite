#pragma once
#include <vector>
#include <string>
#include <Eigen/Core>
#include <set>
#include <utility>
#include "topology.hpp"

int assign_atoms(const std::string& process_atoms,
                 const char* atomname_optional,
                 const std::vector<std::string>& Anames,
                 const std::vector<std::string>& Bnames,
                 const Eigen::MatrixXd& distmat,
                 std::vector<int>& assignBofA,
                 std::vector<int>& assignAofB,
                 double threshold);

int assign_atoms_resinfo(const std::vector<std::string>& Anames,
                         const std::vector<std::string>& Bnames,
                         const std::vector<std::string>& Aresnames,
                         const std::vector<std::string>& Bresnames,
                         const std::vector<int>& Aresids,
                         const std::vector<int>& Bresids,
                         std::vector<int>& assignBofA,
                         std::vector<int>& assignAofB);

void assign_atoms_connectivity(const Eigen::MatrixXd& distmat,
                               const topology& Atop,
                               const topology& Btop,
                               const std::set<std::pair<std::string, std::string> > &forbid_assign,
                               std::vector<int>& assignBofA,
                               std::vector<int>& assignAofB,
                               std::vector<int>& Adepth,
                               double threshold,
                               bool must_be_identical_names);

void unassign_atoms_forbidding(const std::vector<std::string>& Anames,
                               const std::vector<std::string>& Bnames,
                               const std::vector<std::string>& Aresnames,
                               const std::vector<std::string>& Bresnames,
                               const std::set<std::pair<std::string, std::string> > &forbid_assign,
                               std::vector<int>& assignBofA,
                               std::vector<int>& assignAofB);


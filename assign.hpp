#pragma once
#include <vector>
#include <string>
#include <Eigen/Core>
#include "topology.hpp"

int assign_atoms(const std::string& process_atoms,
                 const char* atomname_optional,
                 const std::vector<std::string>& Anames,
                 const std::vector<std::string>& Bnames,
                 const Eigen::MatrixXd& distmat,
                 std::vector<int>& assignBofA,
                 std::vector<int>& assignAofB,
                 double threshold);

void assign_atoms_connectivity(const Eigen::MatrixXd& distmat,
                               const topology& Atop,
                               const topology& Btop,
                               std::vector<int>& assignBofA,
                               std::vector<int>& assignAofB,
                               std::vector<int>& Adepth,
                               double threshold,
                               bool must_be_identical_names);

#pragma once
#include <vector>
#include <string>
#include <Eigen/Core>
#include <set>
#include <utility>
#include "topology.hpp"

enum assigner_action {
    ASSIGN_ACCEPT,
    ASSIGN_REJECT,
    ASSIGN_NR
};

struct assigner_conditions {
    assigner_action action;
    std::string res1, atom1, res2, atom2;
};

typedef std::vector<assigner_conditions> assigner_dictionary;

int assign_atoms(const std::string& process_atoms,
                 const char* atomname_optional,
                 const topology &Atop,
                 const topology &Btop,
                 const Eigen::MatrixXd& distmat,
                 std::vector<int>& assignBofA,
                 std::vector<int>& assignAofB,
                 double threshold);

int assign_atoms_resinfo(const topology &Atop,
                         const topology &Btop,
                         const assigner_dictionary& dictionary,
                         std::vector<int>* assignBofA,
                         std::vector<int>* assignAofB);

void assign_atoms_connectivity(const topology& Atop,
                               const topology& Btop,
                               const Eigen::MatrixXd& distmat,
                               std::vector<int>& assignBofA,
                               std::vector<int>& assignAofB,
                               std::vector<int>& Adepth,
                               double threshold,
                               bool must_be_identical_names);


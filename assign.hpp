#pragma once
#include <vector>
#include <string>
#include <Eigen/Core>

void assign_atoms(const std::string& process_atoms,
                  const std::vector<std::string>& Anames,
                  const std::vector<std::string>& Bnames,
                  const Eigen::MatrixXd& distmat,
                  std::vector<int>& assignBofA,
                  std::vector<int>& assignAofB,
                  double threshold);

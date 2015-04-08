#pragma once

#include <vector>
#include <string>
#include <Eigen/Core>

Eigen::VectorXd set_selected_mass(const std::vector<std::string>& atomnames,
                  const std::string& sel);
void fit_selected(const Eigen::VectorXd& Amass,
                  const Eigen::VectorXd& Bmass,
                  Eigen::Matrix3Xd& Acoords,
                  Eigen::Matrix3Xd& Bcoords);


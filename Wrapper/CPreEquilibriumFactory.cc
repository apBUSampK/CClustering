/**
 * Copyright (c) 2026 Savva Savenkov, Artemii Novikov
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy of this software and
 * associated documentation files (the "Software"), to deal in the Software without restriction,
 * including without limitation the rights to use, copy, modify, merge, publish, distribute,
 * sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all copies or
 * substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT
 * NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
 * DAMAGES OR OTHER* OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#include "CPreEquilibriumFactory.hh"

#include "CoordinateMSTClustering.hh"

#include <COLA.hh>

#include <memory>
#include <stdexcept>
#include <string>
#include <unordered_map>

using namespace cola;

std::unique_ptr<VFilter> CPreEquilibriumFactory::Create(const std::unordered_map<std::string, std::string>& param_map) {
  if (param_map.at("clustering_type") == "GMST") {
    if (auto it = param_map.find("stat_exen_type"); it != param_map.end()) {
      excitation_energy_type_ = std::stoi(it->second);
    }
    if (auto it = param_map.find("consider_coulomb"); it != param_map.end()) {
      repulsion_ = std::stoi(it->second);
    }
    if (auto it = param_map.find("simulate_momentum"); it != param_map.end()) {
      momentum_ = std::stoi(it->second);
    }
    return std::make_unique<CoordinateMSTClustering>(repulsion_.value_or(false), momentum_.value_or(true),
                                                     excitation_energy_type_.value_or(4));
  }
  throw std::runtime_error("Clustering type is unrecognized");
}
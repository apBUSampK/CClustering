//
// Created by apbus_amp_k on 26.02.24.
//

#include <memory>

#include "CPreEquilibriumFactory.hh"
#include "MSTClustering.hh"
#include "CPreEquilibrium.hh"

cola::VFilter* CPreEquiilibriumFactory::create(const std::map<std::string, std::string>& paramMap) {
  if (paramMap.at("clustering_type") == "GMST") {
    if (auto it = paramMap.find("stat_exen_type"); it != paramMap.end())
    {
      excitationEnergyType = std::stoi(it->second);
    }
    if (auto it = paramMap.find("consider_coulomb"); it != paramMap.end())
    {
      repulsion = std::stoi(it->second);
    }
    if (auto it = paramMap.find("simulate_momentum"); it != paramMap.end())
    {
      momentum = std::stoi(it->second);
    }
    return new CoordinateMSTClustering(repulsion.value_or(false), momentum.value_or(true), excitationEnergyType.value_or(7));
  } else {
    throw std::runtime_error("Clustering type is unrecognized");
  }
}
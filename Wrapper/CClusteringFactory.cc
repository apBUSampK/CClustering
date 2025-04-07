//
// Created by apbus_amp_k on 26.02.24.
//

#include <memory>

#include "CClusteringFactory.hh"
#include "MSTClustering.h"

cola::VFilter* CClusteringFactory::create(const std::map<std::string, std::string>& paramMap) {
  if (paramMap.at("clustering_type") == "GMST") {
    uint stat_exen_type = std::stoi(paramMap.at("stat_exen_type"));
    uint consider_rep = std::stoi(paramMap.at("consider_coulomb"));
    return new GMSTClustering(stat_exen_type, consider_rep);
  } else {
    throw std::runtime_error("Clustering type is unrecognized");
  }
}
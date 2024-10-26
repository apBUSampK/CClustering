//
// Created by apbus_amp_k on 26.02.24.
//

#include <memory>

#include "CClusteringFactory.hh"
#include "MSTClustering.h"

cola::VFilter* CClusteringFactory::create(const std::map<std::string, std::string>& paramMap) {
  // extract neede parameters
  // if (paramMap.at("clustering_type") == "GMST") {
  //   return new GMSTClustering(paramMap.at("readerID"));
  // } else {
  //   throw std::runtime_error("Clustering type is unrecognized");
  // }
}
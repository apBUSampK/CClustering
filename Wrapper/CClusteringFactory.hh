//
// Created by apbus_amp_k on 26.02.24.
//

#ifndef CCLUSTERING_CCLUSTERINGFACTORY_HH
#define CCLUSTERING_CCLUSTERINGFACTORY_HH

#include "COLA.hh"

class CClusteringFactory final: public cola::VFactory {
public:
    cola::VFilter* create(const std::map<std::string, std::string>&) final;

};

#endif //CGLAUBER_CCLUSTERINGFACTORY_HH

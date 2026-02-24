//
// Created by apbus_amp_k on 26.02.24.
//

#ifndef CCLUSTERING_CCLUSTERINGFACTORY_HH
#define CCLUSTERING_CCLUSTERINGFACTORY_HH

#include "COLA.hh"

#include <optional>

class CPreEquiilibriumFactory final: public cola::VFactory {
public:
    cola::VFilter* create(const std::map<std::string, std::string>&) final;
private:
    std::optional<bool> repulsion;
    std::optional<bool> momentum;
    std::optional<int> excitationEnergyType;
};

#endif //CGLAUBER_CCLUSTERINGFACTORY_HH

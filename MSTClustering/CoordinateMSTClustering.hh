#ifndef CCLUSTERING_GMSTCLUSTERING_H
#define CCLUSTERING_GMSTCLUSTERING_H

#include <memory>

#include <algorithm>

#include "MSTClustering.hh"
#include "G4ExcitationHandler.hh"
#include "G4FermiPhaseSpaceDecay.hh"
#include "G4ReactionProductVector.hh"
#include "G4NucleiProperties.hh"

class CoordinateMSTClustering : public MSTClustering {
  public:
    CoordinateMSTClustering() = delete;
    CoordinateMSTClustering(bool consider_rep, bool extra_momentum, int stat_exen_type) : _consider_rep(consider_rep),
      _extra_momentum(extra_momentum), _stat_exen_type(stat_exen_type) {};
  private:

    static constexpr double nucleonAverMass = 0.93891875434*CLHEP::GeV;
    bool _consider_rep;
    bool _extra_momentum;
    uint _stat_exen_type;

    std::vector<Edge> get_edges(const cola::EventData&) final;
    std::unique_ptr<cola::EventData> get_clusters(std::unique_ptr<cola::EventData>&&) final;
    
    cola::EventParticles calculate_momentum(std::vector<std::vector<cola::Particle*>> noMomClusters, double ExEnA, double ExEnB, CLHEP::Hep3Vector boostA, CLHEP::Hep3Vector boostB, cola::EventParticles rnucsA, cola::EventParticles rnucsB, std::vector<int> rmapsA, std::vector<int> rmapsB);
    cola::EventParticles _process_side(const cola::EventData&, cola::ParticleClass);

    G4FermiPhaseSpaceDecay phaseSpaceDecay;
};

// helper functions
cola::LorentzVector ToColaLorentzVector(G4LorentzVector& lv);

#endif //CCLUSTERING_GMSTCLUSTERING_H
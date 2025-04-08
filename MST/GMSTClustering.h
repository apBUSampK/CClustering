#ifndef CCLUSTERING_GMSTCLUSTERING_H
#define CCLUSTERING_GMSTCLUSTERING_H

#include <memory>
#include <algorithm>

#include "MSTClustering.h"
#include "G4ExcitationHandler.hh"
#include "G4SystemOfUnits.hh"
#include "G4FermiPhaseSpaceDecay.hh"
#include "G4ReactionProductVector.hh"
#include "G4NucleiProperties.hh"

class GMSTClustering : public MSTClustering {
  public:
    GMSTClustering() = default;
    GMSTClustering(uint stat_exen_type, uint consider_rep) : stat_exen_type_(stat_exen_type), consider_rep_(static_cast<bool>(consider_rep)) {};
  protected:
    using edge = std::pair<double, std::pair<uint, uint>>;
  private:

    static constexpr double nucleonAverMass = 0.93891875434*CLHEP::GeV;
    uint sourceA;
    uint sourceAb;
    uint A = 0;
    uint Ab = 0;
    uint Z = 0;
    uint Zb = 0;
    uint stat_exen_type_ = 4;
    bool consider_rep_ = true;
    cola::EventParticles pls_;

    std::vector<edge> get_vertices(const cola::EventData&) final;
    std::unique_ptr<cola::EventData> get_clusters(std::unique_ptr<cola::EventData>&&, const Node&) final;
    double get_cd(double Ex, uint A);
    std::pair<double, double> get_exens();
    std::vector<std::vector<uint>> get_comps(double, cola::ParticleClass);
    void dfs(std::shared_ptr<Node> node, std::vector<bool>& visited, std::vector<uint>& component, double cd, cola::ParticleClass pClass);
    std::vector<cola::Particle*> fragments_from_clusters(const std::vector<std::vector<uint>>&);
    cola::EventParticles calculate_momentum(std::vector<std::vector<cola::Particle*>> noMomClusters, double ExEnA, double ExEnB, CLHEP::Hep3Vector boostA, CLHEP::Hep3Vector boostB, cola::EventParticles rnucsA, cola::EventParticles rnucsB, std::vector<int> rmapsA, std::vector<int> rmapsB);
    CLHEP::Hep3Vector get_boost(CLHEP::Hep3Vector p, double E, double A);
    cola::LorentzVector ToColaLorentzVector(G4LorentzVector lv);


	G4double eps0 = 2.17 * MeV;
	G4double alphaPow = -1.02;
	G4double d0 = 2.7;
	G4double aColRel = 5.315;
	G4double aSurfRel  = 0.017;
	G4double aVol     = 0.054;

    G4FermiPhaseSpaceDecay phaseSpaceDecay;

    G4double SpecAa = 0;
    G4double SpecAb = 0;

    G4double a_opt = 2.243;
    G4double b_opt = 3.183 * MeV;
    G4double c_opt = 0.99;
    G4double d_opt = 0.29041;
};

#endif //CCLUSTERING_GMSTCLUSTERING_H
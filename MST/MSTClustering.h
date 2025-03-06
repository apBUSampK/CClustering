//
// Created by _amp_ on 10/21/24.
//

#ifndef CCLUSTERING_MSTCLUSTERING_H
#define CCLUSTERING_MSTCLUSTERING_H

#include <memory>
#include <algorithm>

#include "COLA.hh"

#include "G4ExcitationHandler.hh"
#include "G4SystemOfUnits.hh"
#include "G4FermiPhaseSpaceDecay.hh"
#include "G4ReactionProductVector.hh"
#include "G4NucleiProperties.hh"


class MSTClustering : public cola::VConverter {
public:
    MSTClustering() = default;
    ~MSTClustering() override = default;
    MSTClustering& operator=(const MSTClustering&) = delete; //verbosity, may delete later
    MSTClustering& operator=(MSTClustering&&) = delete;
    MSTClustering(const MSTClustering&) = delete;
    MSTClustering(MSTClustering&&) = delete;

    std::unique_ptr<cola::EventData> operator()(std::unique_ptr<cola::EventData>&& data) final {
        construct_tree(get_vertices(*data), data->particles.size());
        return get_clusters(std::move(data), *tree.at(0));
    }

protected:
    using iPair = std::pair<uint, uint>;        // vertices pair
    using edge = std::pair<double, iPair>;      // edge

    // a single Node with children.
    class Node {
    public:
        explicit Node(uint vertice) : height(0.), vertices(1, vertice), children() {}
        Node() : Node(0) {}
        Node(std::shared_ptr<Node>&& first, std::shared_ptr<Node>&& second, double height) : height(height),
        vertices(first->vertices), children(std::make_pair(first, second)) {
            vertices.insert(vertices.end(), second->vertices.begin(), second->vertices.end()); // append second vector
        }

        [[nodiscard]] double get_height() const {return height;}
        [[nodiscard]] std::vector<uint> get_vertices() const { return vertices;}
        [[nodiscard]] std::pair<std::shared_ptr<Node>, std::shared_ptr<Node>> get_children() const {return children;}

    private:
        double height;
        std::vector<uint> vertices;
        std::pair<std::shared_ptr<Node>, std::shared_ptr<Node>> children;

    };

    std::vector<std::shared_ptr<Node>> tree;

private:

    void construct_tree(std::vector<edge>&& verticeData, size_t size);

    // get full graph data from EventData
    virtual std::vector<edge> get_vertices(const cola::EventData&) = 0;
    // construct clusters (Node is MST root)
    virtual std::unique_ptr<cola::EventData> get_clusters(std::unique_ptr<cola::EventData>&&, const Node&) = 0;
};

class GMSTClustering : public MSTClustering {
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

    std::vector<edge> get_vertices(const cola::EventData&) final;
    std::unique_ptr<cola::EventData> get_clusters(std::unique_ptr<cola::EventData>&&, const Node&) final;
    double get_cd(double Ex, uint A);
    std::pair<double, double> get_exens(uint statType);
    std::vector<std::vector<uint>> get_connected_components(double);
    void dfs(std::shared_ptr<Node> node, std::vector<bool>& visited, std::vector<uint>& component, double cd);
    std::vector<cola::Particle*> fragments_from_clusters(const std::vector<std::vector<uint>>&, const cola::EventParticles&);
    cola::EventParticles calculate_momentum(std::vector<std::vector<cola::Particle*>> noMomClusters, double ExEnA, double ExEnB, CLHEP::Hep3Vector boostA, CLHEP::Hep3Vector boostB, cola::EventParticles rnucsA, cola::EventParticles rnucsB, std::vector<int> rmapsA, std::vector<int> rmapsB);
    CLHEP::Hep3Vector get_boost(uint pZ, uint A);
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

#endif //CCLUSTERING_MSTCLUSTERING_H

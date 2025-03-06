#ifndef REPULSION_HH
#define REPULSION_HH

#include <cmath>
#include <limits>
#include <iostream>
#include <vector>
#include <memory>

#include "COLA.hh"
#include "G4SystemOfUnits.hh"
#include "G4ReactionProductVector.hh"

namespace RepulsionStage {

constexpr double fm = 1e-15 * CLHEP::m;
constexpr double theta = 0.3;
constexpr double totalTime = 200 * fm / CLHEP::c_light;
constexpr double iterations = 1000;
constexpr double max_adaptive_delta = std::numeric_limits<double>::max();

std::vector<cola::Particle*> CalculateRepulsion(std::vector<cola::Particle*> frags, cola::EventParticles nucleons, const std::vector<int>& maps);

class BHNode {
 public:
  int totalA; // total proton count
  G4ThreeVector cr; // mean coordinates of the charges in box
  G4ThreeVector ctr; // coordinates of the box center
  std::vector<std::unique_ptr<BHNode>> children; // child nodes
  int index; // -1 if > 1 particles, index in nucleons vector otherwise
  double size; // size of the box

  BHNode() = default;
  BHNode(double size, G4ThreeVector ctr) : size(size), ctr(ctr), totalA(0), cr({0.0, 0.0, 0.0}), index(-1) {};
  ~BHNode() = default;
  void Divide();
};

class BHTree {
 public:
  explicit BHTree(const cola::EventParticles* nucleons, std::vector<cola::Particle*>* frags, const std::vector<int>* maps);

  std::vector<G4ThreeVector> Iterate(double time_delta);

  double GetAdaptiveTimeDelta() const;

 private:
  std::unique_ptr<BHNode> rootnode_;
  std::vector<cola::Particle*>* frags_;
  const std::vector<int>* maps_;
  std::vector<G4ThreeVector> fs_;

  void BuildBHTree(const cola::EventParticles* nucleons);

  std::unique_ptr<BHNode> InitializeRoot(const cola::EventParticles* nucleons);

  void GetForces(const BHNode* node);

  G4ThreeVector Force(const BHNode* rootnode, const BHNode* node) const;

  G4ThreeVector DuoForce(const G4ThreeVector vec, const double& from_totalA) const;

  void InsertNucleon(const std::unique_ptr<BHNode>& node, const G4ThreeVector& cords, int pIndex);
};

}

#endif
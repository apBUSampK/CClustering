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

cola::EventParticles CalculateRepulsion(cola::EventParticles&& frags);

class BHNode {
 public:
  int Z; // total charge
  cola::Vector3<double> cr; // mean coordinates of the charges in box
  cola::Vector3<double> ctr; // coordinates of the box center
  std::vector<std::unique_ptr<BHNode>> children; // child nodes
  int index; // -1 if > 1 particles, index in nucleons vector otherwise
  double size; // size of the box

  BHNode() = default;
  BHNode(double size, cola::Vector3<double>& ctr) : size(size), ctr(ctr), Z(0), cr({0.0, 0.0, 0.0}), index(-1) {};
  ~BHNode() = default;
  void Divide();
};

class BHTree {
 public:
  explicit BHTree(cola::EventParticles& frags);

  std::vector<cola::Vector3<double>> Iterate(double time_delta);

  double GetAdaptiveTimeDelta() const;

 private:
  std::unique_ptr<BHNode> rootnode_;
  cola::EventParticles& frags_;
  std::vector<cola::Vector3<double>> fs_;

  void BuildBHTree(const cola::EventParticles& frags);

  std::unique_ptr<BHNode> InitializeRoot(const cola::EventParticles& frags);

  void GetForces(const BHNode* node);

  cola::Vector3<double> Force(const BHNode* rootnode, const BHNode* node) const;

  cola::Vector3<double> DuoForce(const cola::Vector3<double> vec, const double& from_totalA) const;

  void InsertNucleon(const std::unique_ptr<BHNode>& node, const cola::Vector3<double>& cords, int pIndex);
};

}

#endif
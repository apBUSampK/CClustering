/**
 * Copyright (c) 2026 Savva Savenkov, Artemii Novikov
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy of this software and
 * associated documentation files (the "Software"), to deal in the Software without restriction,
 * including without limitation the rights to use, copy, modify, merge, publish, distribute,
 * sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all copies or
 * substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT
 * NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
 * DAMAGES OR OTHER* OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#ifndef REPULSION_HH
#define REPULSION_HH

#include <CLHEP/Units/PhysicalConstants.h>
#include <CLHEP/Units/SystemOfUnits.h>
#include <EventData.hh>
#include <LorentzVector.hh>

#include <cmath>
#include <limits>
#include <memory>
#include <vector>

namespace repulsion_stage {

  constexpr double kFm = 1e-15 * CLHEP::m;
  constexpr double kTheta = 0.3;
  constexpr double kTotalTime = 200 * kFm / CLHEP::c_light;
  constexpr double kIterations = 1000;
  constexpr double kMaxAdaptiveDelta = std::numeric_limits<double>::max();

  cola::EventParticles CalculateRepulsion(cola::EventParticles&& frags);

  class BHNode {
   public:
    int z;                                          // total charge
    int n_part;                                     // number of particles in the node
    cola::Vector3<double> cr;                       // mean coordinates of the charges in box
    cola::Vector3<double> ctr;                      // coordinates of the box center
    std::vector<std::unique_ptr<BHNode>> children;  // child nodes
    int index;                                      // -1 if > 1 particles, index in nucleons vector otherwise
    double size;                                    // size of the box

    BHNode() = default;
    BHNode(double size, const cola::Vector3<double>& ctr)
        : size(size), ctr(ctr), z(0), n_part(0), cr({.x = 0.0, .y = 0.0, .z = 0.0}), index(-1) {};
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
    std::vector<cola::Vector3<double>> fs_;  // Forces

    void BuildBHTree(const cola::EventParticles& frags);

    static std::unique_ptr<BHNode> InitializeRoot(const cola::EventParticles& frags);

    void GetForces(const BHNode* node);

    cola::Vector3<double> Force(const BHNode* rootnode, const BHNode* node) const;

    static cola::Vector3<double> DuoForce(cola::Vector3<double> vec, double from_z);

    void InsertFragment(const std::unique_ptr<BHNode>& node, const cola::Vector3<double>& cords, int p_index, int z);
  };

}  // namespace repulsion_stage

#endif
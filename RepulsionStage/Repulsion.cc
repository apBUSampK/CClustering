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

#include "Repulsion.hh"

#include <CLHEP/Units/PhysicalConstants.h>
#include <EventData.hh>
#include <LorentzVector.hh>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <memory>
#include <vector>

using namespace cola;

namespace repulsion_stage {

  cola::EventParticles CalculateRepulsion(cola::EventParticles&& frags) {
    if (frags.empty()) {
      return frags;
    }

    double time = 0.0;
    double delta_time = kTotalTime / static_cast<double>(kIterations);

    while (time < kTotalTime) {
      BHTree bhtree(frags);

      double temp_timedelta = std::min(delta_time, bhtree.GetAdaptiveTimeDelta());

      auto r_delta = bhtree.Iterate(temp_timedelta);

      for (size_t i = 0; i < frags.size(); i++) {
        frags[i].position.t += temp_timedelta;
        frags[i].position.x += r_delta[i].x;
        frags[i].position.y += r_delta[i].y;
        frags[i].position.z += r_delta[i].z;
      }
      time += temp_timedelta;
    }
    return frags;
  }

  BHTree::BHTree(cola::EventParticles& frags) : frags_(frags), fs_(frags.size(), {.0, .0, .0}) {
    BuildBHTree(frags);
    GetForces(rootnode_.get());
  }

  std::unique_ptr<BHNode> BHTree::InitializeRoot(const cola::EventParticles& frags) {
    cola::Vector3 min_vector;
    cola::Vector3 max_vector;
    min_vector.x = std::ranges::min_element(frags, [](const cola::Particle& left, const cola::Particle& right) {
                     return left.position.x < right.position.x;
                   })->position.x;
    min_vector.y = std::ranges::min_element(frags, [](const cola::Particle& left, const cola::Particle& right) {
                     return left.position.y < right.position.y;
                   })->position.y;
    min_vector.z = std::ranges::min_element(frags, [](const cola::Particle& left, const cola::Particle& right) {
                     return left.position.z < right.position.z;
                   })->position.z;
    max_vector.x = std::ranges::max_element(frags, [](const cola::Particle& left, const cola::Particle& right) {
                     return left.position.x < right.position.x;
                   })->position.x;
    max_vector.y = std::ranges::max_element(frags, [](const cola::Particle& left, const cola::Particle& right) {
                     return left.position.y < right.position.y;
                   })->position.y;
    max_vector.z = std::ranges::max_element(frags, [](const cola::Particle& left, const cola::Particle& right) {
                     return left.position.z < right.position.z;
                   })->position.z;

    auto range = max_vector - min_vector;
    double max_range = std::max({range.x, range.y, range.z});
    return std::make_unique<BHNode>(max_range, (min_vector + max_vector) / 2);
  }

  void BHTree::BuildBHTree(const cola::EventParticles& frags) {
    rootnode_ = InitializeRoot(frags);

    for (size_t i = 0; i < frags.size(); i++) {
      InsertFragment(rootnode_, frags.at(i).position.SpatialPart(), i, frags.at(i).GetAZ().second);
    }
  }

  void BHTree::InsertFragment(const std::unique_ptr<BHNode>& node, const cola::Vector3<double>& cords, int p_index,
                              int fragment_z) {
    if (node->n_part == 0) {
      node->total_z += fragment_z;
      node->cr = cords;
      node->index = p_index;
      node->n_part++;
      return;
    }
    if (node->n_part == 1) {
      node->Divide();

      int i = 0;
      if (node->cr.x > node->ctr.x) {
        i |= 1;
      }
      if (node->cr.y > node->ctr.y) {
        i |= 2;
      }
      if (node->cr.z > node->ctr.z) {
        i |= 4;
      }
      InsertFragment(node->children[i], node->cr, node->index, node->total_z);

      i = 0;
      if (cords.x > node->ctr.x) {
        i |= 1;
      }
      if (cords.y > node->ctr.y) {
        i |= 2;
      }
      if (cords.z > node->ctr.z) {
        i |= 4;
      }
      InsertFragment(node->children[i], cords, p_index, fragment_z);

      node->cr = (node->cr + cords) / 2;
      node->total_z += fragment_z;
      node->index = -1;
      node->n_part++;
      return;
    }
    node->cr = (node->cr * node->total_z + cords) / (node->total_z + fragment_z);
    node->total_z += fragment_z;
    node->n_part++;

    int i = 0;
    if (cords.x > node->ctr.x) {
      i |= 1;
    }
    if (cords.y > node->ctr.y) {
      i |= 2;
    }
    if (cords.z > node->ctr.z) {
      i |= 4;
    }
    InsertFragment(node->children[i], cords, p_index, fragment_z);
  }

  double BHTree::GetAdaptiveTimeDelta() const {
    double min_time = kMaxAdaptiveDelta;
    for (size_t i = 0; i < frags_.size(); i++) {
      auto mnt = frags_.at(i).momentum;
      double mval = std::sqrt(mnt.x * mnt.x + mnt.y * mnt.y + mnt.z * mnt.z);
      if ((frags_.at(i).GetAZ().second == 0u) || (fs_[i].Mag() == 0.0) || (mval == 0.0)) {
        continue;
      }
      min_time = std::min(min_time, 0.05 * mval / fs_[i].Mag());
    }
    return min_time;
  }

  std::vector<cola::Vector3<double>> BHTree::Iterate(double time_delta) {
    std::vector<cola::Vector3<double>> r_delta(frags_.size());

    for (size_t i = 0; i < frags_.size(); ++i) {
      if (frags_.at(i).GetAZ().second == 0) {
        continue;
      }
      cola::LorentzVector momentum = frags_.at(i).momentum;
      cola::Vector3<double> half_dp = fs_[i] * time_delta * 0.5;

      cola::Vector3<double> mid_v =
          (momentum.SpatialPart() + half_dp) / std::sqrt((momentum.SpatialPart() + half_dp).Mag2() + momentum.Mag2());
      r_delta[i] =
          cola::Vector3<double>{.x = time_delta * mid_v.x, .y = time_delta * mid_v.y, .z = time_delta * mid_v.z};

      cola::Vector3<double> pvec = momentum.SpatialPart() + 2 * half_dp;
      momentum.x = pvec.x;
      momentum.y = pvec.y;
      momentum.z = pvec.z;
      momentum.e = std::sqrt((momentum.SpatialPart() + 2 * half_dp).Mag2() + momentum.Mag2());

      frags_.at(i).momentum = momentum;
    }
    return r_delta;
  }

  void BHTree::GetForces(const BHNode* node) {
    if (node->total_z == 0) {
      return;
    }
    if (node->n_part == 1) {
      fs_[node->index] += Force(rootnode_.get(), node);
      return;
    }
    for (const auto& child : node->children) {
      GetForces(child.get());
    }
  }

  cola::Vector3<double> BHTree::Force(const BHNode* rootnode, const BHNode* node) const {
    if ((rootnode->total_z == 0) || ((rootnode->index != -1) && rootnode->index == node->index)) {
      return {.x = 0.0, .y = 0.0, .z = 0.0};
    }
    if (rootnode->n_part == 1) {
      return DuoForce(node->cr - rootnode->cr, rootnode->total_z);
    }
    if ((rootnode->size / (node->cr - rootnode->cr).Mag()) < kTheta) {
      return DuoForce(node->cr - rootnode->cr, rootnode->total_z);
    }
    cola::Vector3<double> total_force = {.x = 0.0, .y = 0.0, .z = 0.0};
    for (const auto& child : rootnode->children) {
      total_force += Force(child.get(), node);
    }
    return total_force;
  }

  cola::Vector3<double> BHTree::DuoForce(const cola::Vector3<double> vec, const double from_z) {
    return vec * (CLHEP::elm_coupling * from_z / std::pow(vec.Mag(), 3));
  }

  void BHNode::Divide() {
    children.resize(8);
    for (size_t i = 0; i < 8; i++) {
      cola::Vector3<double> offset{.x = (((i & 1) != 0u) ? 1 : -1) * size / 4.0,
                                   .y = (((i & 2) != 0u) ? 1 : -1) * size / 4.0,
                                   .z = (((i & 4) != 0u) ? 1 : -1) * size / 4.0};
      children[i] = std::make_unique<BHNode>(size / 2.0, ctr + offset);
    }
  }
}  // namespace repulsion_stage
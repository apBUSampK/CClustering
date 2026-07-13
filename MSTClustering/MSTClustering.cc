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

#include "MSTClustering.hh"

#include "ExcitationEnergy.hh"
#include "Repulsion.hh"

#include <G4ExcitationHandler.hh>
#include <G4NucleiProperties.hh>

#include <limits>
#include <stack>

using namespace cola;

void MSTClustering::construct_trees(std::vector<Edge>&& edgeData, std::vector<Node>& nodes) {
  std::unordered_map<cola::Particle*, Node*> treeA;
  std::unordered_map<cola::Particle*, Node*> treeB;
  treeA.reserve(std::distance(spectIterA, spectIterB));
  treeB.reserve(std::distance(spectIterB, endIter));

  // no resizes are allowed (otherwise segfault)
  nodes.reserve(2 * (treeA.bucket_count() + treeB.bucket_count()));

  // initialize trees
  for (auto iter = spectIterA; iter < spectIterB; ++iter) {
    nodes.emplace_back(*iter);
    treeA.emplace(&(*iter), &nodes.back());
  }
  for (auto iter = spectIterB; iter < endIter; ++iter) {
    nodes.emplace_back(*iter);
    treeB.emplace(&(*iter), &nodes.back());
  }

  // sort the edges for hierarchical tree
  std::sort(edgeData.begin(), edgeData.end(), [](const Edge& l, const Edge& r) { return l.size < r.size; });

  // merge nodes into complete trees using edgeData
  for (const auto& edge : edgeData) {
    auto v1 = edge.vert.first;
    auto v2 = edge.vert.second;
    switch (edge.p_class) {
      case cola::ParticleClass::kSpectatorA:
        if (treeA[v1] != treeA[v2]) {
          nodes.emplace_back(treeA[v1], treeA[v2], edge.size);
          for (const auto vertex : nodes.back().vertices) {
            treeA[vertex] = &nodes.back();
          }
        }
        break;
      case cola::ParticleClass::kSpectatorB:
        if (treeB[v1] != treeB[v2]) {
          nodes.emplace_back(treeB[v1], treeB[v2], edge.size);
          for (const auto vertex : nodes.back().vertices) {
            treeB[vertex] = &nodes.back();
          }
        }
        break;
      default:
        throw(std::logic_error("An edge between non-Spectators have been formed!"));
    }
  }

  // set root nodes

  rootA = !treeA.empty() ? treeA.begin()->second : nullptr;
  rootB = !treeB.empty() ? treeB.begin()->second : nullptr;
}

std::unique_ptr<cola::EventData> MSTClustering::operator()(std::unique_ptr<cola::EventData>&& data) {
  // sort particles by p_class (spectators are last)
  std::sort(data->particles.begin(), data->particles.end(),
            [](const cola::Particle& l, const cola::Particle& r) { return l.p_class < r.p_class; });
  spectIterA = std::find_if(data->particles.begin(), data->particles.end(),
                            [](cola::Particle p) { return p.p_class == cola::ParticleClass::kSpectatorA; });
  spectIterB = std::find_if(data->particles.begin(), data->particles.end(),
                            [](cola::Particle p) { return p.p_class == cola::ParticleClass::kSpectatorB; });
  endIter = data->particles.end();
  // construct trees
  std::vector<Node> nodes;  // stores all nodes in trees (improves cpu cache hits)
  construct_trees(get_edges(*data), nodes);
  // divide trees and process resulting pre-fragments
  return get_clusters(std::move(data));
}

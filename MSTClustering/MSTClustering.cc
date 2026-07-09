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

#include <EventData.hh>

#include <algorithm>
#include <iterator>
#include <memory>
#include <stdexcept>
#include <unordered_map>
#include <utility>
#include <vector>

using namespace cola;

void MSTClustering::construct_trees(std::vector<Edge>&& edgeData, std::vector<Node>& nodes) {
  std::unordered_map<cola::Particle*, Node*> tree_a;
  std::unordered_map<cola::Particle*, Node*> tree_b;
  tree_a.reserve(std::distance(spectIterA, spectIterB));
  tree_b.reserve(std::distance(spectIterB, endIter));

  // no resizes are allowed (otherwise segfault)
  nodes.reserve(2 * (tree_a.bucket_count() + tree_b.bucket_count()));

  // initialize trees
  for (auto iter = spectIterA; iter < spectIterB; ++iter) {
    nodes.emplace_back(*iter);
    tree_a.emplace(&(*iter), &nodes.back());
  }
  for (auto iter = spectIterB; iter < endIter; ++iter) {
    nodes.emplace_back(*iter);
    tree_b.emplace(&(*iter), &nodes.back());
  }

  // sort the edges for hierarchical tree
  std::ranges::sort(edgeData, [](const Edge& l, const Edge& r) { return l.size < r.size; });

  // merge nodes into complete trees using edgeData
  for (const auto& edge : edgeData) {
    auto* v1 = edge.vert.first;
    auto* v2 = edge.vert.second;
    switch (edge.p_class) {
      case cola::ParticleClass::kSpectatorA:
        if (tree_a[v1] != tree_a[v2]) {
          nodes.emplace_back(tree_a[v1], tree_a[v2], edge.size);
          for (auto* const vertex : nodes.back().vertices) {
            tree_a[vertex] = &nodes.back();
          }
        }
        break;
      case cola::ParticleClass::kSpectatorB:
        if (tree_b[v1] != tree_b[v2]) {
          nodes.emplace_back(tree_b[v1], tree_b[v2], edge.size);
          for (auto* const vertex : nodes.back().vertices) {
            tree_b[vertex] = &nodes.back();
          }
        }
        break;
      default:
        throw(std::logic_error("An edge between non-Spectators have been formed!"));
    }
  }

  // set root nodes

  rootA = !tree_a.empty() ? tree_a.begin()->second : nullptr;
  rootB = !tree_b.empty() ? tree_b.begin()->second : nullptr;
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

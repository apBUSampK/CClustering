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

#ifndef CCLUSTERING_MSTCLUSTERING_H
#define CCLUSTERING_MSTCLUSTERING_H

#include <COLA.hh>
#include <EventData.hh>

#include <memory>
#include <optional>
#include <utility>
#include <vector>

namespace cola {

  class MSTClustering : public VConverter {
   public:
    MSTClustering() = default;
    ~MSTClustering() override = default;
    MSTClustering& operator=(const MSTClustering&) = delete;
    MSTClustering& operator=(MSTClustering&&) = delete;
    MSTClustering(const MSTClustering&) = delete;
    MSTClustering(MSTClustering&&) = delete;

    std::unique_ptr<EventData> operator()(std::unique_ptr<EventData>&& data) final;

   protected:
    using nPair = std::pair<Particle*, Particle*>;  // vertices pair
    struct Edge {                                   // edge
      nPair vert;
      double size;
      ParticleClass p_class;
      Edge(nPair vert, double size, ParticleClass p_class) : vert(std::move(vert)), size(size), p_class(p_class) {};
    };

    // a single Node with children.
    struct Node {
      explicit Node(Particle& vertex) : vertices({&vertex}) {}
      Node() = default;
      Node(Node* first, Node* second, double height)
          : height(height), vertices(first->vertices), children(std::make_pair(first, second)) {
        vertices.insert(vertices.end(), second->vertices.begin(), second->vertices.end());  // append second vector
      }

      Node(Node&&) noexcept = default;

      double height = 0.;
      std::vector<Particle*> vertices;
      std::optional<std::pair<Node*, Node*>> children;
    };

    Node* root_a_;
    Node* root_b_;

    // it is reasonable to use the iterators for quick access to spectators even at this abstract level, so this class
    // has it
    EventParticles::iterator spect_iter_a_;
    EventParticles::iterator spect_iter_b_;
    EventParticles::iterator end_iter_;

   private:
    void ConstructTrees(std::vector<Edge>&& edge_data, std::vector<Node>& nodes);

    // get full graph data from EventData
    virtual std::vector<Edge> GetEdges(const EventData&) = 0;
    // construct clusters
    virtual std::unique_ptr<EventData> GetClusters(std::unique_ptr<EventData>&&) = 0;
  };
}  // namespace cola
#endif  // CCLUSTERING_MSTCLUSTERING_H

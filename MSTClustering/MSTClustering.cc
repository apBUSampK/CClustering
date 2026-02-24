//
// Created by _amp_ on 10/21/24.
//

#include "MSTClustering.hh"
#include "G4NucleiProperties.hh"
#include "Repulsion.hh"
#include "G4ExcitationHandler.hh"
#include "ExcitationEnergy.hh"
#include <limits>
#include <stack>

void MSTClustering::construct_trees(std::vector<Edge> &&edgeData)
{
    std::map<cola::Particle*, std::shared_ptr<Node>> treeA;
    std::map<cola::Particle*, std::shared_ptr<Node>> treeB;

    // initialize trees
    for (auto iter = spectIterA; iter < spectIterB; iter++)
        treeA.emplace(&(*iter), std::make_shared<Node>(&(*iter)));
    for (auto iter = spectIterB; iter < endIter; iter++)
        treeB.emplace(&(*iter), std::make_shared<Node>(&(*iter)));
        
    // merge nodes into complete trees using edgeData
    for (auto it : edgeData)
    {
        auto v1 = it.vert.first;
        auto v2 = it.vert.second;
        switch (it.pClass)
        {
        case cola::ParticleClass::spectatorA:
            if (treeA[v1].get() != treeA[v2].get()) {
                auto newNode = std::make_shared<Node>(std::move(treeA[v1]), std::move(treeA[v2]), it.size);
                for (const auto vertice : newNode->vertices) treeA[vertice] = newNode;
            }
            break;
        case cola::ParticleClass::spectatorB:
            if (treeB[v1].get() != treeB[v2].get()) {
                auto newNode = std::make_shared<Node>(std::move(treeB[v1]), std::move(treeB[v2]), it.size);
                for (const auto vertice : newNode->vertices) treeB[vertice] = newNode;
            }
            break;
        default:
            throw(std::logic_error("An edge between non-Spectators have been formed!"));
        }
    }

    // set root nodes

    rootA = treeA.begin()->second;
    rootB = treeB.begin()->second;
}

//
// Created by _amp_ on 10/21/24.
//

#include "MSTClustering.h"


void MSTClustering::construct_tree(std::vector<edge>&& verticeData, size_t size) {
    tree.clear();
    for (int i = 0; i < size; i++)
        tree.push_back(std::make_shared<Node>(i));
    std::sort(verticeData.begin(), verticeData.end());
    for (auto it: verticeData) {
        uint v1 = it.second.first;
        uint v2 = it.second.second;
        if (tree[v1].get() != tree[v2].get()) {
            auto newNode = std::make_shared<Node>(std::move(tree[v1]), std::move(tree[v2]), it.first);
            for (const auto vertice: newNode->get_vertices())
                tree[vertice] = newNode;
        }
    }
}
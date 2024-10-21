//
// Created by _amp_ on 10/21/24.
//

#ifndef CCLUSTERING_MSTCLUSTERING_H
#define CCLUSTERING_MSTCLUSTERING_H

#include <memory>
#include <algorithm>

#include "COLA.hh"


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
    using iPair = std::pair<uint, uint>;
    using edge = std::pair<double, iPair>;

private:
    class Node {
    public:
        explicit Node(uint vertice) : height(0.), vertices(1, vertice), children() {}
        Node() : Node(0) {}
        Node(std::shared_ptr<Node>&& first, std::shared_ptr<Node>&& second, double height) : height(height),
        vertices(first->vertices), children(std::make_pair(first, second)) {
            vertices.insert(vertices.end(), second->vertices.begin(), second->vertices.end()); // append second vector
        }

    private:
        double height;
        std::vector<uint> vertices;
        std::pair<std::shared_ptr<Node>, std::shared_ptr<Node>> children;

    };

    std::vector<std::shared_ptr<Node>> tree;

    void construct_tree(std::vector<edge>&& verticeData, size_t size) {
        tree.clear();
        for (int i = 0; i < size; i++)
            tree.push_back(std::make_shared<Node>(i));
        std::sort(verticeData.begin(), verticeData.end());
        for (auto it: verticeData) {
            uint v1 = it.second.first;
            uint v2 = it.second.second;
            if (tree[v1].get() != tree[v2].get()) {
                tree[v1] = std::make_shared<Node>(std::move(tree[v1]), std::move(tree[v2]), it.first);
                tree[v2] = tree[v1];
            }
        }
    }

    virtual std::vector<edge> get_vertices(const cola::EventData&) = 0;
    virtual std::unique_ptr<cola::EventData> get_clusters(std::unique_ptr<cola::EventData>&&, const Node&) = 0;
};


#endif //CCLUSTERING_MSTCLUSTERING_H

#pragma once

#include <vector>
#include <cassert>

struct Edge {
    uint32_t source;
    uint32_t target;

    Edge(uint32_t s, uint32_t t)
        : source(s)
        , target(t)
    { }
};

struct Graph {
    size_t num_vertices;
    std::vector<Edge> edges;
    std::vector<uint32_t> nOffsets;
    std::vector<uint32_t> neighbours;

public: /* Methods: */

    Graph(size_t n, std::vector<Edge> es);

    uint32_t numberOfNeighbours(uint32_t i) const {
        assert (i < num_vertices);
        assert (nOffsets[i] <= nOffsets[i+1]);
        return nOffsets[i + 1] - nOffsets[i];
    }

public: /* Static methods: */

    static Graph repeat(Graph gr, uint32_t k);
    
    static Graph triangles(uint32_t num);
    static Graph lines(uint32_t num);
    static Graph line(uint32_t num);
    static Graph points(uint32_t num);

    static Graph grid(uint32_t width, uint32_t height);
    static Graph hypercube(uint32_t dim);
    static Graph tree(uint32_t levels, uint32_t k = 2);
    static Graph kneser(uint32_t n, uint32_t k);
    static Graph tesselation(uint32_t n);
    static Graph complete(uint32_t n) { return kneser(n, 1); }
    static Graph circle(uint32_t n);
    static Graph ladder(uint32_t height, uint32_t n);

    template <typename Sampler>
    static Graph erModel(uint32_t n, Sampler sampler);

    template <typename UniformSampler>
    static Graph ba1Model(uint32_t n, UniformSampler sample01);
};
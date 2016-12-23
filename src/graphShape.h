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

// TODO: better name
struct GraphShape {
    size_t num_vertices;
    std::vector<Edge> edges;
    std::vector<uint32_t> nOffsets;
    std::vector<uint32_t> neighbours;

public: /* Methods: */

    GraphShape(size_t n, std::vector<Edge> es);

    uint32_t numberOfNeighbours(uint32_t i) const {
        assert (i < num_vertices);
        assert (nOffsets[i] <= nOffsets[i+1]);
        return nOffsets[i + 1] - nOffsets[i];
    }

public: /* Static methods: */

    static GraphShape repeat(GraphShape shape, uint32_t k);
    
    static GraphShape triangles(uint32_t num);
    static GraphShape lines(uint32_t num);
    static GraphShape line(uint32_t num);
    static GraphShape points(uint32_t num);

    static GraphShape grid(uint32_t width, uint32_t height);
    static GraphShape hypercube(uint32_t dim);
    static GraphShape tree(uint32_t levels, uint32_t k = 2);
    static GraphShape kneser(uint32_t n, uint32_t k);
    static GraphShape tesselation(uint32_t n);
    static GraphShape complete(uint32_t n) { return kneser(n, 1); }
    static GraphShape circle(uint32_t n);
    static GraphShape ladder(uint32_t height, uint32_t n);

    template <typename Sampler>
    static GraphShape erModel(uint32_t n, Sampler sampler);

    template <typename UniformSampler>
    static GraphShape ba1Model(uint32_t n, UniformSampler sample01);
};
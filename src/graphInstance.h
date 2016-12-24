#pragma once

#include <vector>

#include "graph.h"
#include "vec2.h"

// TODO: better name
struct GraphInstance {
    const Graph* gr;
    std::vector<vec2f> positions;
    std::vector<float> energies;
    float totalEnergy;

public: /* Methods: */

    GraphInstance(const Graph& s, std::vector<vec2f> verts);

    // TODO: result either delta or new energy
    void update(uint32_t v, vec2f coord);
    static GraphInstance randomised(const Graph& gr);

    float vertSqrDist(uint32_t s, uint32_t t) const;

    float edgeEnergy(Edge e) const;

    float vertexEnergy(uint32_t v) const;

    GraphInstance neighbour(float radius) const;
};
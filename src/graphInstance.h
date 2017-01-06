#pragma once

#include "graph.h"
#include "vec2.h"

#include <vector>



// TODO: better name
struct GraphInstance {
public: /* Types: */

    struct Delta {
        float energyDelta;
        uint32_t vertex;
        vec2f position;
    };

public: /* Fields: */

    const Graph* gr;
    std::vector<vec2f> positions;
    std::vector<float> energies;
    float totalEnergy;

public: /* Methods: */

    GraphInstance(const Graph& s, std::vector<vec2f> verts);
    static GraphInstance randomised(const Graph& gr);

    Delta neighbour(float radius);
    void commit(const Delta& delta);

private: /* Methods: */

    float vertSqrDist(uint32_t s, uint32_t t) const;
    float edgeEnergy(uint32_t source, uint32_t target) const;
    float vertexEnergy(uint32_t v) const;
    float computeEnergyDelta(uint32_t v, vec2f coord);
};
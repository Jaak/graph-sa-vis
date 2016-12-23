#pragma once

#include <vector>
#include "graphShape.h"
#include "vec2.h"

// TODO: better name
struct GraphInstance {
    const GraphShape* shape;
    std::vector<vec2f> positions;
    std::vector<float> energies;
    float totalEnergy;

    GraphInstance(const GraphShape& s, std::vector<vec2f> verts);

    // TODO: result either delta or new energy
    void update(uint32_t v, vec2f coord);
    static GraphInstance randomised(const GraphShape& shape);

    float vertSqrDist(uint32_t s, uint32_t t) const;

    float edgeEnergy(Edge e) const;

    float vertexEnergy(uint32_t v) const;
};


// TODO: move to a member
GraphInstance neighbour(GraphInstance inst, float radius);

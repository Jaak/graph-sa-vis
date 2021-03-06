#include "graph.h"
#include "graphInstance.h"
#include "options.h"
#include "random.h"

#include <iterator>
#include <cassert>
#include <numeric>

template <typename Iter>
typename std::iterator_traits<Iter>::value_type variance(Iter begin, Iter end) {
    const auto size = std::distance(begin, end);
    const auto mean = std::accumulate(begin, end, 0.0f) / size;
    const auto sqrMean = std::inner_product(begin, end, begin, 0.0f) / size;
    const auto variance = sqrMean - mean*mean;
    return variance;
}

GraphInstance::GraphInstance(const Graph& graph, std::vector<vec2f> pos)
    : gr(&graph)
    , positions(std::move(pos))
    , energies(gr->num_vertices, 0.0f)
    , totalEnergy(0)
{
    for (uint32_t i = 0; i < gr->num_vertices; ++ i) {
        const auto e = vertexEnergy(i);
        energies[i] = e;
        totalEnergy += e;
    }
}

float GraphInstance::vertSqrDist(uint32_t s, uint32_t t) const {
    return sqrdist(positions[s], positions[t]);
}

float GraphInstance::edgeEnergy(uint32_t source, uint32_t target) const {
    float energy = 0;

    const auto edgeWeight = EdgeDistanceWeight * gr->num_vertices / gr->edges.size();
    energy += edgeWeight*vertSqrDist(source, target);

    for (size_t j = 0; j < gr->edges.size(); ++ j) {
        const auto e = gr->edges[j];
        if (source == e.source || source == e.target ||
            target == e.source || target == e.target) {
            continue;
        }


        const auto p1 = positions[source];
        const auto q1 = positions[target];
        const auto p2 = positions[e.source];
        const auto q2 = positions[e.target];
        if (intersects(p1, q1, p2, q2)) {
            energy += 0.5f*CrossingWeight;
        }
    }

    return energy;
}

float GraphInstance::vertexEnergy(uint32_t v) const {
    float energy = 0.0f;

    if (ClosenessWeight > 0.0f) {
        const auto vertWeight = 0.5f * ClosenessWeight / gr->num_vertices;
        for (uint32_t i = 0; i < gr->num_vertices; ++ i) {
            if (i != v) {
                energy += vertWeight / (0.001f + vertSqrDist(v, i));
            }
        }
    }

    if (EdgeFromVertexDistanceWeight > 0.0f) {
        for (auto edge : gr->edges) {
            if (edge.source == v || edge.target == v)
                continue;

            const auto sqrDist = sqrDistanceFromLine(
                positions[edge.source], positions[edge.target], positions[v]);
            if (sqrDist) {
                energy += EdgeFromVertexDistanceWeight / (0.001f + *sqrDist);
            }
        }
    }

    const auto begin = gr->nOffsets[v];
    const auto end = gr->nOffsets[v + 1];
    assert(begin <= end);
    const auto n = end - begin;

    for (uint32_t i = begin; i < end; ++ i) {
        energy += 0.5f*edgeEnergy(v, gr->neighbours[i]);
    }

    if (AngleWeight > 0.0f && n >= 2) {
        constexpr float DOUBLE_PI = 2*M_PI;
        const auto p0 = positions[v];
        const auto p1 = positions[gr->neighbours[begin]];
        const auto u = normalised(p1 - p0);

        std::vector<float> angles;
        angles.reserve(n);
        angles.push_back(DOUBLE_PI);
        for (uint32_t j = begin + 1; j < end; ++ j) {
            const auto pj = positions[gr->neighbours[j]];
            const auto v = normalised(pj - p0);
            auto a = std::atan2(cross(u, v) , dot(u, v));
            if (a < 0.0f) a += DOUBLE_PI;
            angles.push_back(a);
        }
        
        std::sort(angles.begin(), angles.end());

        for (uint32_t i = n; i --> 1; ) {
            angles[i] -= angles[i - 1];
        }

        energy += AngleWeight * variance(angles.begin(), angles.end());
    }

    return energy;
}

GraphInstance GraphInstance::randomised(const Graph& gr) {
    std::vector<vec2f> coords;
    coords.reserve(gr.num_vertices);
    for (size_t i = 0; i < gr.num_vertices; ++ i) {
        coords.emplace_back(randFloat(BoxPadding, BoundingBoxWidth - BoxPadding),
                            randFloat(BoxPadding, BoundingBoxHeight - BoxPadding));
    }

    return GraphInstance {gr, coords};
}

// TODO: avoid mutating the state
float GraphInstance::computeEnergyDelta(uint32_t v, vec2f pos) {
    const auto oldPos = positions[v];
    const auto oldEnergy = energies[v];

    positions[v] = pos;
    const auto newEnergy = vertexEnergy(v);
    positions[v] = oldPos;
    const auto delta = newEnergy - oldEnergy;
    return delta;
}

void GraphInstance::commit(const GraphInstance::Delta& delta) {
    positions[delta.vertex] = delta.position;
    energies[delta.vertex] += delta.energyDelta;
    totalEnergy += delta.energyDelta;
}

GraphInstance::Delta GraphInstance::neighbour(float radius) {
    const uint32_t i = randInt(0, positions.size() - 1);
    const auto pos = positions[i];
    for (;;) {
        const auto d = radius * sampleCircle();
        const auto p = pos + d;
        if (p.x >= BoxPadding && p.x + BoxPadding < BoundingBoxWidth &&
            p.y >= BoxPadding && p.y + BoxPadding < BoundingBoxHeight)
        {
            const auto energyDelta = computeEnergyDelta(i, p);
            return { energyDelta, i, p };
        }
    }
}
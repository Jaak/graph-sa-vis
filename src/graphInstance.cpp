#include "options.h"
#include "graphInstance.h"
#include "graphShape.h"
#include "random.h"

#include <iterator>
#include <cassert>
#include <numeric>

template <typename Iter>
typename std::iterator_traits<Iter>::value_type variance(Iter begin, Iter end) {
    const auto size = std::distance(begin, end);
    const auto iSize = 1.0 / size;
    const auto mean = iSize * std::accumulate(begin, end, 0.0);
    const auto sqrMean = iSize * std::inner_product(begin, end, begin, 0.0);
    const auto variance = sqrMean - mean*mean;
    return variance;
}

GraphInstance::GraphInstance(const GraphShape& s, std::vector<vec2f> pos)
    : shape(&s)
    , positions(std::move(pos))
    , energies(shape->num_vertices, 0.0f)
    , totalEnergy(0)
{
    for (uint32_t i = 0; i < shape->num_vertices; ++ i) {
        const auto e = vertexEnergy(i);
        energies[i] = e;
        totalEnergy += e;
    }
}

void GraphInstance::update(uint32_t v, vec2f coord) {
    positions[v] = coord;
    const auto oldEnergy = energies[v];
    const auto newEnergy = vertexEnergy(v);
    energies[v] = newEnergy;
    totalEnergy -= oldEnergy;
    totalEnergy += newEnergy;
}

float GraphInstance::vertSqrDist(uint32_t s, uint32_t t) const {
    return sqrdist(positions[s], positions[t]);
}

float GraphInstance::edgeEnergy(Edge e) const {
    float energy = 0;

    const auto edgeWeight = EdgeDistanceWeight * shape->num_vertices / shape->edges.size();
    energy += edgeWeight*vertSqrDist(e.source, e.target);

    for (size_t j = 0; j < shape->edges.size(); ++ j) {
        const auto e2 = shape->edges[j];
        if (e.source == e2.source || e.source == e2.target ||
            e.target == e2.source || e.target == e2.target) {
            continue;
        }


        const auto p1 = positions[e.source];
        const auto q1 = positions[e.target];
        const auto p2 = positions[e2.source];
        const auto q2 = positions[e2.target];
        if (intersects(p1, q1, p2, q2)) {
            energy += 0.5*CrossingWeight;
        }
    }

    return energy;
}

float GraphInstance::vertexEnergy(uint32_t v) const {
    float energy = 0;

    if (ClosenessWeight > 0) {
        const auto vertWeight = 0.5 * ClosenessWeight / shape->num_vertices;
        for (uint32_t i = 0; i < shape->num_vertices; ++ i) {
            if (i != v) {
                energy += vertWeight / (0.001 + vertSqrDist(v, i));
            }
        }
    }

    if (EdgeFromVertexDistanceWeight > 0) {
        for (auto edge : shape->edges) {
            if (edge.source == v || edge.target == v)
                continue;

            const auto sqrDist = sqrDistanceFromLine(
                positions[edge.source], positions[edge.target], positions[v]);
            if (sqrDist) {
                energy += EdgeFromVertexDistanceWeight / (0.001 + *sqrDist);
            }
        }
    }

    const auto begin = shape->nOffsets[v];
    const auto end = shape->nOffsets[v + 1];
    assert(begin <= end);
    const auto n = end - begin;

    for (uint32_t i = begin; i < end; ++ i) {
        energy += 0.5*edgeEnergy(Edge(v, shape->neighbours[i]));
    }

    if (AngleWeight > 0 && n >= 2) {
        const auto p0 = positions[v];
        const auto p1 = positions[shape->neighbours[begin]];
        const auto u = normalised(p1 - p0);

        std::vector<float> angles;
        angles.reserve(n);
        angles.push_back(2*M_PI);
        for (uint32_t j = begin + 1; j < end; ++ j) {
            const auto pj = positions[shape->neighbours[j]];
            const auto v = normalised(pj - p0);
            auto a = std::atan2(cross(u, v) , dot(u, v));
            if (a < 0) a += 2*M_PI;
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

GraphInstance GraphInstance::randomised(const GraphShape& shape) {
    std::vector<vec2f> coords;
    coords.reserve(shape.num_vertices);
    for (size_t i = 0; i < shape.num_vertices; ++ i) {
        coords.emplace_back(randFloat(BoxPadding, BoundingBoxWidth - BoxPadding),
                            randFloat(BoxPadding, BoundingBoxHeight - BoxPadding));
    }

    return GraphInstance {shape, coords};
}

GraphInstance neighbour(GraphInstance inst, float radius) {
    const auto i = randInt(0, inst.positions.size() - 1);
    auto& vertex = inst.positions[i];
    for (;;) {
        const auto d = radius * sampleCircle();
        const auto p = vertex + d;
        if (p.x >= BoxPadding && p.x + BoxPadding < BoundingBoxWidth &&
            p.y >= BoxPadding && p.y + BoxPadding < BoundingBoxHeight)
        {
            inst.update(i, p);
            break;
        }
    }

    return inst;
}

#include "graphShape.h"


GraphShape::GraphShape(size_t n, std::vector<Edge> es)
    : num_vertices(n)
    , edges(std::move(es))
{
    nOffsets.resize(n + 1, 0);
    neighbours.resize(2 * edges.size(), 0);

    // Count the number of outgoing edges:
    for (auto edge : edges) {
        nOffsets[1 + edge.source] ++;
        nOffsets[1 + edge.target] ++;
    }

    for (uint32_t i = 0; i < n; ++ i) {
        nOffsets[i+1] += nOffsets[i];
    }

    std::vector<uint32_t> index = nOffsets;
    for (auto edge : edges) {
        neighbours[index[edge.source] ++] = edge.target;
        neighbours[index[edge.target] ++] = edge.source;
    }
}

inline std::vector<uint32_t> choose(uint32_t n, uint32_t k) {
    if (k == 0) return {0};
    if (n == 0) return {};
    const auto& xs = choose(n - 1, k);
    const auto& ys = choose(n - 1, k - 1);
    std::vector<uint32_t> result;
    result.reserve(xs.size() + ys.size());
    for (auto x : xs) result.push_back(x << 1);
    for (auto y : ys) result.push_back((y << 1) ^ 1);
    return result;
}

inline std::vector<Edge> mkHypercubeEdges(uint32_t dim) {
    if (dim == 0) {
        return std::vector<Edge>{};
    }

    std::vector<Edge> result;
    const auto es = mkHypercubeEdges(dim - 1);
    const auto n = 1u << (dim - 1);
    for (auto e : es) {
        result.push_back(e);
        result.emplace_back(e.source + n, e.target + n);
    }

    for (uint32_t i = 0; i < n; ++ i) {
        result.emplace_back(i, i + n);
    }

    return result;
}

GraphShape GraphShape::hypercube(uint32_t dim) {
    return GraphShape {1u << dim, mkHypercubeEdges(dim) };
}

GraphShape GraphShape::tesselation(uint32_t n) {
    std::vector<Edge> edges;
    uint32_t prevNode = 0;
    uint32_t currentNode = 1;
    for (uint32_t l = 1; l < n; ++ l) {
        for (uint32_t i = 0; i < l; ++ i) {
            edges.emplace_back(prevNode, currentNode);
            edges.emplace_back(prevNode, currentNode + 1);
            edges.emplace_back(currentNode, currentNode + 1);
            prevNode ++;
            currentNode ++;
        }

        currentNode ++;
    }

    return GraphShape {currentNode, edges};
}

GraphShape GraphShape::repeat(GraphShape shape, uint32_t k) {
    const uint32_t n = shape.num_vertices;
    const uint32_t num_vertices = n * k;
    std::vector<Edge> edges = std::move(shape.edges);
    const uint32_t ne = edges.size();
    edges.reserve(shape.edges.size() * k);
    for (uint32_t i = 1; i < k; ++ i) {
        for (uint32_t j = 0; j < ne; ++ j) {
            edges.emplace_back(
                edges[j].source + i*n,
                edges[j].target + i*n);
        }
    }

    return GraphShape {num_vertices, edges};
}

GraphShape GraphShape::kneser(uint32_t n, uint32_t k) {
    assert (n > k);

    const auto vec = choose(n, k);
    const uint32_t num_vertices = vec.size();
    std::vector<Edge> edges;
    for (uint32_t i = 0; i < num_vertices; ++ i)
        for (uint32_t j = i + 1; j < num_vertices; ++ j)
            if ((vec[i] & vec[j]) == 0)
                edges.emplace_back(i, j);

    return GraphShape {num_vertices, edges};
}

GraphShape GraphShape::triangles(uint32_t num) {
    return GraphShape::repeat(GraphShape::circle(3), num);
}

GraphShape GraphShape::lines(uint32_t num) {
    const uint32_t num_vertices = 2 * num;
    std::vector<Edge> edges;
    for (uint32_t i = 0; i < num_vertices; i += 2) {
        edges.emplace_back(i, i + 1);
    }

    return GraphShape {num_vertices, edges};
}

GraphShape GraphShape::line(uint32_t num) {
    const uint32_t num_vertices = num;
    std::vector<Edge> edges;
    for (uint32_t i = 1; i < num_vertices; ++ i) {
        edges.emplace_back(i - 1, i);
    }

    return GraphShape {num_vertices, edges};
}

GraphShape GraphShape::points(uint32_t num) {
    return GraphShape {num, std::vector<Edge>{}};
}

GraphShape GraphShape::tree(uint32_t levels, uint32_t k) {
    uint32_t prev = 0;
    uint32_t current = 1;
    uint32_t levelSize = 1;
    std::vector<Edge> edges;
    for (uint32_t l = 0; l < levels; ++ l) {
        for (uint32_t j = 0; j < levelSize; ++ j) {
            for (size_t i = 0; i < k; ++ i)
                edges.emplace_back(prev, current ++);

            prev ++;
        }

        levelSize = current - prev;
    }
    
    return GraphShape {current, edges};
}

GraphShape GraphShape::grid(uint32_t width, uint32_t height) {
    assert (width > 0 && height > 0);

    const auto mkCoord = [width](uint32_t x, uint32_t y) {
        return x + width * y;
    };
    
    const uint32_t num_vertices = width * height;
    std::vector<Edge> edges;

    for (uint32_t i = 0; i < width; ++ i) {
        for (uint32_t j = 0; j < height; ++ j) {
            if (i > 0) edges.emplace_back(mkCoord(i, j), mkCoord(i - 1, j));
            if (j > 0) edges.emplace_back(mkCoord(i, j), mkCoord(i, j - 1));
        }
    }

    return GraphShape {num_vertices, edges};
}

GraphShape GraphShape::circle(uint32_t n) {
    std::vector<Edge> edges;
    for (uint32_t i = 0; i < n; ++ i)
        edges.emplace_back(i, (i + 1) % n);
    
    return GraphShape {n, edges};
}

GraphShape GraphShape::ladder(uint32_t height, uint32_t n) {
    const uint32_t num_vertices = height * n;
    std::vector<Edge> edges;
    for (uint32_t h = 0; h < height; ++ h) {
        for (uint32_t i = 0; i < n; ++ i) {
            edges.emplace_back(i + h*n, (i + 1) % n + h*n);
            if (h > 0)
                edges.emplace_back(i + (h - 1)*n, i + h*n);
        }
    }

    return GraphShape {num_vertices, edges};
}

// Barabási–Albert model (with m = 1 and m0 = 1)
template <typename UniformSampler>
GraphShape GraphShape::ba1Model(uint32_t n, UniformSampler sample01) {

    std::vector<uint32_t> degrees (n, 0);
    std::vector<Edge> edges;
    uint32_t sumDegrees = 0;

    // Add first two nodes:
    edges.emplace_back(0, 1);
    sumDegrees = 2;
    degrees[0] = 1;
    degrees[1] = 1;

    // Add rest of the nodes:
    for (uint32_t i = 2; i < n; ++ i) {
        const auto x = sample01();
        uint32_t acc = 0;
        for (uint32_t j = 0; j < i; ++ j) {
            acc += degrees[j];
            if (x * sumDegrees <= acc) {
                edges.emplace_back(j, i);
                sumDegrees += 2;
                degrees[i] += 1;
                degrees[j] += 1;
                break;
            }
        }
    }

    return GraphShape {n, edges};
}

// Erdős–Rényi model
template <typename Sampler>
GraphShape GraphShape::erModel(uint32_t n, Sampler sampler) {
    std::vector<Edge> edges;
    for (uint32_t i = 0; i < n; ++ i) {
        for (uint32_t j = i + 1; j < n; ++ j) {
            if (sampler()) {
                edges.emplace_back(i, j);
            }
        }
    }

    return GraphShape {n, edges};
}
#pragma once

#include <vector>

struct Edge {
    uint32_t source;
    uint32_t target;

    Edge(uint32_t s, uint32_t t)
        : source(s)
        , target(t)
    { }
};

struct GraphShape {
    size_t num_vertices;
    std::vector<Edge> edges;
    std::vector<uint32_t> nOffsets;
    std::vector<uint32_t> neighbours;

public: /* Methods: */

    GraphShape(size_t n, std::vector<Edge> es)
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
};

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

inline GraphShape GraphShape::hypercube(uint32_t dim) {
    return GraphShape {1u << dim, mkHypercubeEdges(dim) };
}

inline GraphShape GraphShape::tesselation(uint32_t n) {
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

inline GraphShape GraphShape::repeat(GraphShape shape, uint32_t k) {
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

inline GraphShape GraphShape::kneser(uint32_t n, uint32_t k) {
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

inline GraphShape GraphShape::triangles(uint32_t num) {
    return GraphShape::repeat(GraphShape::circle(3), num);
}

inline GraphShape GraphShape::lines(uint32_t num) {
    const uint32_t num_vertices = 2 * num;
    std::vector<Edge> edges;
    for (uint32_t i = 0; i < num_vertices; i += 2) {
        edges.emplace_back(i, i + 1);
    }

    return GraphShape {num_vertices, edges};
}

inline GraphShape GraphShape::line(uint32_t num) {
    const uint32_t num_vertices = num;
    std::vector<Edge> edges;
    for (uint32_t i = 1; i < num_vertices; ++ i) {
        edges.emplace_back(i - 1, i);
    }

    return GraphShape {num_vertices, edges};
}

inline GraphShape GraphShape::points(uint32_t num) {
    return GraphShape {num, std::vector<Edge>{}};
}

inline GraphShape GraphShape::tree(uint32_t levels, uint32_t k) {
    uint32_t prev = 0;
    uint32_t current = 1;
    uint32_t levelSize = 1;
    std::vector<Edge> edges;
    for (uint32_t l = 0; l < levels; ++ l) {
        for (uint j = 0; j < levelSize; ++ j) {
            for (size_t i = 0; i < k; ++ i)
                edges.emplace_back(prev, current ++);

            prev ++;
        }

        levelSize = current - prev;
    }
    
    return GraphShape {current, edges};
}

inline GraphShape GraphShape::grid(uint32_t width, uint32_t height) {
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

inline GraphShape GraphShape::circle(uint32_t n) {
    std::vector<Edge> edges;
    for (uint32_t i = 0; i < n; ++ i)
        edges.emplace_back(i, (i + 1) % n);
    
    return GraphShape {n, edges};
}

inline GraphShape GraphShape::ladder(uint32_t height, uint32_t n) {
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
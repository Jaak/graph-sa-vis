#include <iostream>
#include <vector>
#include <cassert>
#include <cairo/cairo.h>
#include <random>
#include <tuple>
#include <cmath>
#include <algorithm>
#include <iterator>
#include <thread>
#include <functional>
#include <future>

#include "vec2.h"
#include "graphShape.h"

using namespace std::chrono;

// Distances given in centimeters:
constexpr float BoundingBoxWidth = 30;
constexpr float BoundingBoxHeight = 20;

// SA Weights:
constexpr float CrossingWeight = 0.3;
constexpr float ClosenessWeight = 100;
constexpr float EdgeDistanceWeight = 1;

constexpr float EdgeFromVertexDistanceWeight = 0;
constexpr float AngleWeight = 0;

// SA parameters:
constexpr float GammaFactor = 0.98;
constexpr int StepFactor = 200;
constexpr float InitialTemperature = 1e6;
constexpr float MinimumMoveRadius = BoundingBoxHeight / 200;

/************************
 * Controller interface *
 ************************/

#if 0

class Controller;
class GraphShape;

GraphShape* ctrlMake();
void ctrlDestroy(GraphShape* ctrl);
void ctrlLoadGraphShape(Controller* ctrl, const GraphShape* shape);
void ctrlRandomiseNodes(Controller* ctrl);
void ctrlStartSA(Controller* ctrl);
void ctrlPauseSA(Controller* ctrl);
void ctrlHandleGraphDraw(Controller* ctrl, cairo_t* cr, int width, int height);
void ctrlHandleEnergyPlotDraw(Controller* ctrl, ciaro_t* cr, int width, int height);

#endif

/***************************
 * Random number generator *
 ***************************/

std::random_device rd;
std::mt19937 gen{rd()};

inline float randFloat(float low, float high) {
    std::uniform_real_distribution<float> dist(low, high);
    return dist(gen);
}

inline vec2f sampleCircle() {
    for (;;) {
        const auto x1 = randFloat(-1, 1);
        const auto x2 = randFloat(-1, 1);
        const auto x1Sqr = x1 * x1;
        const auto x2Sqr = x2 * x2;
        const auto sqrSum = x1Sqr + x2Sqr;
        if (sqrSum < 1.0f) {
            const auto invSqrSum = 1.0f / sqrSum;
            return {(x1Sqr - x2Sqr) * invSqrSum, 2*x1*x2*invSqrSum };
        }
    }
}

inline int randInt(int low, int high) {
    std::uniform_int_distribution<int> dist(low, high);
    return dist(gen);
}

/*********************
 * Some math support *
 *********************/

template <typename Iter>
typename std::iterator_traits<Iter>::value_type variance(Iter begin, Iter end) {
    const auto size = std::distance(begin, end);
    const auto iSize = 1.0 / size;
    const auto mean = iSize * std::accumulate(begin, end, 0.0);
    const auto sqrMean = iSize * std::inner_product(begin, end, begin, 0.0);
    const auto variance = sqrMean - mean*mean;
    return variance;
}

/***********************
 * Instance of a graph *
 ***********************/

struct GraphInstance;

struct GraphInstance {
    const GraphShape* shape;
    std::vector<vec2f> vertices;
    std::vector<float> energies;
    float totalEnergy;

    GraphInstance(const GraphShape& s, std::vector<vec2f> verts)
        : shape(&s)
        , vertices(std::move(verts))
        , energies(shape->num_vertices, 0.0f)
        , totalEnergy(0)
    {
        for (uint32_t i = 0; i < shape->num_vertices; ++ i) {
            const auto e = vertexEnergy(i);
            energies[i] = e;
            totalEnergy += e;
        }
    }

    void update(uint32_t v, vec2f coord) {
        vertices[v] = coord;
        const auto oldEnergy = energies[v];
        const auto newEnergy = vertexEnergy(v);
        energies[v] = newEnergy;
        totalEnergy -= oldEnergy;
        totalEnergy += newEnergy;
    }

    static GraphInstance randomised(const GraphShape& shape);

    inline float vertSqrDist(uint32_t s, uint32_t t) const {
        return sqrdist(vertices[s], vertices[t]);
    }

    float edgeEnergy(Edge e) const {
        float energy = 0;

        const auto edgeWeight = EdgeDistanceWeight * shape->num_vertices / shape->edges.size();
        energy += edgeWeight*vertSqrDist(e.source, e.target);

        for (size_t j = 0; j < shape->edges.size(); ++ j) {
            const auto e2 = shape->edges[j];
            if (e.source == e2.source || e.source == e2.target ||
                e.target == e2.source || e.target == e2.target) {
                continue;
            }


            const auto p1 = vertices[e.source];
            const auto q1 = vertices[e.target];
            const auto p2 = vertices[e2.source];
            const auto q2 = vertices[e2.target];
            if (intersects(p1, q1, p2, q2)) {
                energy += 0.5*CrossingWeight;
            }
        }

        return energy;
    }

    float vertexEnergy(uint32_t v) const {
        float energy = 0;

        if (ClosenessWeight > 0) {
            const auto vertWeight = 0.5 * ClosenessWeight / shape->num_vertices;
            for (uint32_t i = 0; i < shape->num_vertices; ++ i) {
                if (i != v) {
                    energy += vertWeight / (0.001 + vertSqrDist(v, i)); // TODO
                }
            }
        }

        if (EdgeFromVertexDistanceWeight > 0) {
            for (auto edge : shape->edges) {
                if (edge.source == v || edge.target == v)
                    continue;

                const auto sqrDist = sqrDistanceFromLine(
                    vertices[edge.source], vertices[edge.target], vertices[v]);
                if (sqrDist) {
                    energy += EdgeFromVertexDistanceWeight / (0.001 + *sqrDist); // TODO
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
            const auto p0 = vertices[v];
            const auto p1 = vertices[shape->neighbours[begin]];
            const auto u = normalised(p1 - p0);

            std::vector<float> angles;
            angles.reserve(n);
            angles.push_back(2*M_PI);
            for (uint32_t j = begin + 1; j < end; ++ j) {
                const auto pj = vertices[shape->neighbours[j]];
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
};


GraphInstance GraphInstance::randomised(const GraphShape& shape) {
    std::uniform_real_distribution<float> xDist(0, BoundingBoxWidth);
    std::uniform_real_distribution<float> yDist(0, BoundingBoxHeight);

    std::vector<vec2f> coords;
    coords.reserve(shape.num_vertices);
    for (size_t i = 0; i < shape.num_vertices; ++ i) {
        coords.emplace_back(xDist(gen), yDist(gen));
    }

    return GraphInstance {shape, coords};
}

float acceptanceP(float e0, float e1, float T) {
    if (e1 < e0)
        return 1;
    else
        return std::exp(-(e1 - e0)/T);
}

GraphInstance neighbour(GraphInstance inst, float radius) {
    const auto i = randInt(0, inst.vertices.size() - 1);
    auto& vertex = inst.vertices[i];
    for (;;) {
        // const auto a = randFloat(0, 2 * M_PI);
        const auto d = radius * sampleCircle(); // vec2f(std::cos(a), std::sin(a));
        const auto p = vertex + d;
        if (p.x >= 0 && p.x < BoundingBoxWidth && p.y >= 0 && p.y < BoundingBoxHeight) {
            inst.update(i, p);
            break;
        }
    }

    return inst;
}

template <typename T>
using event_type = std::pair<time_point<steady_clock>, T>;

using anneal_result =
    std::pair<GraphInstance, std::vector<event_type<float>>> ;

std::vector<event_type<float>> trackEvents(float& e, std::atomic_bool& running) {
    const size_t samplesPerSecond = 4;
    const auto delayInMilliseconds = 100 / samplesPerSecond;

    std::vector<event_type<float>> energies;
    while (running) {
        std::this_thread::sleep_for(milliseconds(delayInMilliseconds));
        const auto energy = e;
        const auto timeNow = steady_clock::now();
        energies.emplace_back(timeNow, energy);
    }

    return energies;
}

anneal_result anneal(GraphInstance s) {
    auto e = s.totalEnergy;
    auto T = InitialTemperature;
    const auto n = s.shape->num_vertices;
    float R = 2.0f * std::min(BoundingBoxHeight, BoundingBoxWidth) / std::sqrt(n);

    std::cout << "InitialEnergy: " << e << std::endl;

    std::atomic_bool running {true};

    std::thread printerThread {
        [&e, &running]() {
            while (running) {
                std::this_thread::sleep_for(milliseconds(50));
                std::cout << "Energy:        " << e << '\r';
                std::cout.flush();
            }
        }
    };

    auto energiesHandle = std::async(std::launch::async, trackEvents, std::ref(e), std::ref(running));

    const auto limit = n * StepFactor;
    size_t stageCount = 0;
    const auto startTime = steady_clock::now();

    for (;;) {
        uint32_t n = 0;
        float mean = 0.0f;
        float M2 = 0.0f;
        for (uint step = 0; step < limit; ++ step) {
            const auto s1 = neighbour(s, R);
            const auto e1 = s1.totalEnergy;
            if (acceptanceP(e, e1, T) >= randFloat(0, 1)) {
                s = s1;
                e = e1;
            }

            n += 1;
            const auto delta = e - mean;
            mean += delta / n;
            const auto delta2 = e - mean;
            M2 += delta*delta2;
        }

        if (M2 / (n - 1) < 0.01f) {
            break;
        }

        T = T * GammaFactor;
        R = std::max(R - 0.5f, MinimumMoveRadius);
        ++ stageCount;
    }

    running = false;
    const auto endTime = steady_clock::now();
    const duration<double> timeDelta = endTime - startTime;

    printerThread.join();
    const auto energies = energiesHandle.get();

    std::cout << "StageCount:    " << stageCount << "        " << std::endl;
    std::cout << "HistorySize:   " << energies.size() << std::endl;
    std::cout << "FinalEnergy:   " << e << std::endl;
    std::cout << "Time:          " << timeDelta.count() << "s" << std::endl;
    return {s, energies};
}

anneal_result anneal(const GraphShape& shape) {
    return anneal(GraphInstance::randomised(shape));
}

/***************
 * Plot energy *
 ***************/

void plotEnergy(const std::vector<event_type<float>>& events,
                int width, int height, std::string filename)
{
    if (events.size() < 3) {
        return;
    }

    const auto startTime = events.front().first;
    const auto endTime = events.back().first;
    const float xPadding = 0.1;
    const float yPadding = 0.1;

    using time_resolution = std::milli;
    const duration<float, time_resolution> timeSpan = endTime - startTime;

    const auto getXCoord = [&](const event_type<float>& event) -> float {
        const duration<float, time_resolution> t = event.first - startTime;
        return (t / timeSpan) * (1.0f - xPadding) + xPadding / 2.0f;
    };

    const auto minmax = std::minmax_element(events.begin(), events.end(),
        [](const event_type<float>& x, const event_type<float>& y){
            return x.second < y.second;
        });
    const auto minEnergy = minmax.first->second;
    const auto maxEnergy = minmax.second->second;

    const auto getYCoord = [&](const event_type<float>& event) -> float {
        const auto t = 1.0f - (event.second - minEnergy) / (maxEnergy - minEnergy);
        return t * (1.0f - yPadding) + yPadding / 2.0f;
    };

    const auto surface = cairo_image_surface_create(CAIRO_FORMAT_ARGB32, width, height);
    const auto cr = cairo_create(surface);

    cairo_save(cr);
    cairo_scale(cr, width, height);
    
    // Draw white background:
    cairo_set_source_rgb(cr, 1, 1, 1);
    cairo_paint(cr);

    // Draw plot:
    for (const auto& event : events) {
        const auto x = getXCoord(event);
        const auto y = getYCoord(event);
        cairo_line_to(cr, x, y);
    }

    cairo_set_source_rgb(cr, 0, 0, 0);
    cairo_set_line_width(cr, 1.0 / width);
    cairo_stroke(cr);
    cairo_restore(cr);

    // Draw to image:
    cairo_surface_write_to_png(surface, filename.c_str());
    cairo_destroy(cr);
    cairo_surface_destroy(surface);
}

/*****************
 * Graph drawing *
 *****************/

void drawGraphInstance(const GraphInstance& inst,
                       int width, int height, std::string filename)
{
    const auto surface = cairo_image_surface_create(CAIRO_FORMAT_ARGB32, width, height);
    const auto cr = cairo_create(surface);

    cairo_save(cr);
    cairo_scale(cr, width / BoundingBoxWidth, height / BoundingBoxHeight);

    // Draw white background:
    cairo_set_source_rgb(cr, 1, 1, 1);
    cairo_paint(cr);

    // Draw lines:
    cairo_set_source_rgb(cr, 0, 0, 0);
    for (auto edge : inst.shape->edges) {
        const auto p1 = inst.vertices[edge.source];
        const auto p2 = inst.vertices[edge.target];
        cairo_move_to(cr, p1.x, p1.y);
        cairo_line_to(cr, p2.x, p2.y);
    }

    cairo_set_line_width(cr, 0.03);
    cairo_stroke(cr);

    cairo_set_line_cap(cr, CAIRO_LINE_CAP_ROUND);

    // Draw points:
    for (auto vert : inst.vertices) {
        cairo_move_to(cr, vert.x, vert.y);
	    cairo_close_path(cr);
    }

    cairo_set_line_width(cr, 0.2);
    cairo_stroke(cr);
    cairo_restore(cr);

    // Draw to image:
    cairo_surface_write_to_png(surface, filename.c_str());
    cairo_destroy(cr);
    cairo_surface_destroy(surface);
}


int main() {
    const auto shape = GraphShape::grid(10, 10);
    // const auto shape = GraphShape::tesselation(10);
    // const auto shape = GraphShape::tree(5, 2);
    // const auto shape = GraphShape::kneser(5, 2);
    // const auto shape = GraphShape::circle(6);

    // const auto shape = GraphShape::lines(100);
    // const auto shape = GraphShape::points(200);
    // const auto shape = GraphShape::repeat(GraphShape::line(5), 200);
    // const auto shape = GraphShape::repeat(GraphShape::circle(6), 200);
    // const auto shape = GraphShape::ladder(2, 7); 

    const auto result = anneal(shape);
    drawGraphInstance(result.first, 800, 600, "out.png");
    plotEnergy(result.second, 800, 400, "plot.png");
    return 0;
}
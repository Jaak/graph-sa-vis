#pragma once

#include "graph.h"
#include "graphInstance.h"

#include <chrono>

template <typename T>
using event_type = std::pair<
    std::chrono::time_point<std::chrono::steady_clock>, T>;

using anneal_result =
    std::pair<GraphInstance, std::vector<event_type<float>>> ;

anneal_result anneal(GraphInstance s);
anneal_result anneal(const Graph& gr);
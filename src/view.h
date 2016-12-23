#pragma once

#include "annealing.h"

#include <cairo/cairo.h>

void plotEnergy(const std::vector<event_type<float>>& events,
                cairo_surface_t* surface, int width, int height);

void drawGraphInstance(const GraphInstance& inst,
                       cairo_surface_t* surface, int width, int height);
#include "options.h"
#include "view.h"

using namespace std::chrono;

void plotEnergy(const std::vector<event_type<float>>& events,
                cairo_surface_t* surface, int width, int height)
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

    const auto cr = cairo_create(surface);

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

    // Draw to image:
    cairo_surface_flush(surface);
    cairo_destroy(cr);
}

void drawGraphInstance(const GraphInstance& inst,
                       cairo_surface_t* surface, int width, int height)
{
    const auto cr = cairo_create(surface);

    cairo_scale(cr, width / BoundingBoxWidth, height / BoundingBoxHeight);

    // Draw white background:
    cairo_set_source_rgb(cr, 1, 1, 1);
    cairo_paint(cr);

    // Draw lines:
    cairo_set_source_rgb(cr, 0, 0, 0);
    for (auto edge : inst.gr->edges) {
        const auto p1 = inst.positions[edge.source];
        const auto p2 = inst.positions[edge.target];
        cairo_move_to(cr, p1.x, p1.y);
        cairo_line_to(cr, p2.x, p2.y);
    }

    cairo_set_line_width(cr, 0.03);
    cairo_stroke(cr);

    cairo_set_line_cap(cr, CAIRO_LINE_CAP_ROUND);

    // Draw points:
    for (auto vert : inst.positions) {
        cairo_move_to(cr, vert.x, vert.y);
	    cairo_close_path(cr);
    }

    cairo_set_line_width(cr, 0.2);
    cairo_stroke(cr);

    // Draw to image:
    cairo_surface_flush(surface);
    cairo_destroy(cr);
}
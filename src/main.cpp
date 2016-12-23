#include <cairo/cairo.h>

#include "graphShape.h"
#include "annealing.h"
#include "view.h"


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


int main() {
    // const auto shape = GraphShape::grid(10, 10);
    // const auto shape = GraphShape::tesselation(10);
    // const auto shape = GraphShape::tree(5, 2);
    // const auto shape = GraphShape::kneser(9, 1);
    // const auto shape = GraphShape::circle(6);
    // const auto shape = GraphShape::hypercube(4);
    
    // const auto shape = GraphShape::erModel(200, []() -> bool {
    //     return randFloat(0.0f, 1.0f) < 0.006f;
    // });

    // const auto shape = GraphShape::ba1Model(100, []() -> float {
    //     return randFloat(0.0f, 1.0f);
    // });

    // const auto shape = GraphShape::lines(100);
    // const auto shape = GraphShape::points(200);
    // const auto shape = GraphShape::repeat(GraphShape::line(5), 200);
    const auto shape = GraphShape::repeat(GraphShape::circle(4), 5);
    // const auto shape = GraphShape::ladder(2, 7); 

    const auto result = anneal(shape);

    {
        const auto width = 800;
        const auto height = 600;
        const auto surface = cairo_image_surface_create(CAIRO_FORMAT_ARGB32, width, height);
        drawGraphInstance(result.first, surface, width, height);
        cairo_surface_write_to_png(surface, "out.png");
        cairo_surface_destroy(surface);
    }

    {
        const auto width = 800;
        const auto height = 400;
        const auto surface = cairo_image_surface_create(CAIRO_FORMAT_ARGB32, width, height);
        plotEnergy(result.second, surface, width, height);
        cairo_surface_write_to_png(surface, "plot.png");
        cairo_surface_destroy(surface);
    }

    return 0;
}
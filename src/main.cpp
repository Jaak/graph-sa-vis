#include "annealing.h"
#include "graph.h"
#include "view.h"

#include <cairo/cairo.h>


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
    // const auto gr = Graph::grid(10, 10);
    // const auto gr = Graph::tesselation(10);
    // const auto gr = Graph::tree(5, 2);
    // const auto gr = Graph::kneser(9, 1);
    // const auto gr = Graph::circle(6);
    // const auto gr = Graph::hypercube(4);
    
    // const auto gr = Graph::erModel(200, []() -> bool {
    //     return randFloat(0.0f, 1.0f) < 0.006f;
    // });

    // const auto gr = Graph::ba1Model(100, []() -> float {
    //     return randFloat(0.0f, 1.0f);
    // });

    // const auto gr = Graph::lines(100);
    // const auto gr = Graph::points(200);
    // const auto gr = Graph::repeat(Graph::line(3), 100);
    const auto gr = Graph::repeat(Graph::circle(3), 20);
    // const auto gr = Graph::ladder(2, 7); 

    const auto result = anneal(gr);

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
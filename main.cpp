#include "./classes/vector.h"
#include "./classes/voronoi.h"
#include "./classes/optimal_transport.h"
#include "./classes/fluid.h"

int main()
{
    Fluid F(100);
    F.runfluid();
    return 0;

    std::vector<Vector> points(300);
    std::vector<double> lambdas(300, 1);
    double sum = 0;
    for (int i = 0; i < points.size(); i++)
    {
        points[i] = Vector(rand()/(double) RAND_MAX, rand()/(double) RAND_MAX);
        lambdas[i] = rand()/((double) RAND_MAX);
        sum += lambdas[i];
    }
    for (int i = 0; i < points.size(); i++)
    {
        lambdas[i] /= sum;
    }

    //Voronoi voronoi(points, lambdas);
    //voronoi.compute();
    //voronoi.save("voronoi.svg");
    OptimalTransport OT(points, lambdas);
    OT.solve();
    OT.solution.save("voronoi.svg");
}
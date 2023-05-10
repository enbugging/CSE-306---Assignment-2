#include "./classes/vector.h"
#include "./classes/voronoi.h"

int main()
{
    std::vector<Vector> points(30);
    for (int i = 0; i < points.size(); i++)
    {
        points[i] = Vector(rand()/(double) RAND_MAX, rand()/(double) RAND_MAX);
    }

    Voronoi voronoi(points);
    voronoi.compute();
    voronoi.save("voronoi.svg");
}
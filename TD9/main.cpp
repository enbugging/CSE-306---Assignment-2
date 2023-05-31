#include "./classes/edge.h"
#include "./classes/vector.h"
#include "./classes/mesh.h"

#include <map>
#include <vector>
#include <cmath>
#include <chrono>
#include <iostream>

int main() {
    // time
    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

    TriangleMesh mesh;
    mesh.readOBJ("./goethe.obj");
    mesh.center_scale();
    std::map<Edge, std::vector<int>> edgeToTriangleIndices;
    std::map<int, std::vector<int>> vertexToTriangleIndices;
    std::vector<bool> is_vertex_on_boundary(mesh.vertices.size(), false);

    for (int  i = 0; i < mesh.indices.size(); i++)
    {
        int vtxi = mesh.indices[i].vtxi;
        int vtxj = mesh.indices[i].vtxj;
        int vtxk = mesh.indices[i].vtxk;

        edgeToTriangleIndices[Edge(vtxi, vtxj)].push_back(i);
        edgeToTriangleIndices[Edge(vtxj, vtxk)].push_back(i);
        edgeToTriangleIndices[Edge(vtxk, vtxi)].push_back(i);

        vertexToTriangleIndices[vtxi].push_back(i);
        vertexToTriangleIndices[vtxj].push_back(i);
        vertexToTriangleIndices[vtxk].push_back(i);
    }

    std::vector<Edge> boundary_edges;
    for (auto it = edgeToTriangleIndices.begin(); it != edgeToTriangleIndices.end(); it++)
    {
        if (it->second.size() == 1)
        {
            boundary_edges.push_back(it->first);
            is_vertex_on_boundary[it->first.a] = true;
            is_vertex_on_boundary[it->first.b] = true;
        }
    }

    std::vector<Edge> ordered_boundary_edges(boundary_edges.size());
    ordered_boundary_edges[0] = boundary_edges[0];
    for (int i = 1; i < boundary_edges.size(); i++)
    {
        for (int j = 0; j < boundary_edges.size(); j++)
        {
            if (boundary_edges[j].a == ordered_boundary_edges[i-1].b)
            {
                ordered_boundary_edges[i] = boundary_edges[j];
                break;
            }
        }
    }

    for (int i = 0; i < ordered_boundary_edges.size(); i++)
    {
        double theta = i*2*M_PI/ordered_boundary_edges.size();
        Vector circle_vertex = Vector(cos(theta)*.5 + .5, sin(theta)*.5 + .5, 0);
        mesh.vertices[ordered_boundary_edges[i].a] = circle_vertex;
    }

    for (int iter = 0; iter < 5000; iter++)
    {
        std::vector<Vector> updated_vertices(mesh.vertices.size());
        #pragma omp parallel for
        for (int i = 0; i < mesh.vertices.size(); i++)
        {
            if (is_vertex_on_boundary[i])
            {
                updated_vertices[i] = mesh.vertices[i];
                continue;
            }

            Vector avg_neighbour(0, 0, 0);
            int number_of_neighbours = 2 * vertexToTriangleIndices[i].size();
            // looping over neighbouring triangles
            for (int j = 0; j < vertexToTriangleIndices[i].size(); j++)
            {
                if (mesh.indices[vertexToTriangleIndices[i][j]].vtxi != i)
                    avg_neighbour += mesh.vertices[mesh.indices[vertexToTriangleIndices[i][j]].vtxi];
                if (mesh.indices[vertexToTriangleIndices[i][j]].vtxj != i)
                    avg_neighbour += mesh.vertices[mesh.indices[vertexToTriangleIndices[i][j]].vtxj];
                if (mesh.indices[vertexToTriangleIndices[i][j]].vtxk != i)
                    avg_neighbour += mesh.vertices[mesh.indices[vertexToTriangleIndices[i][j]].vtxk];
            }

            updated_vertices[i] = avg_neighbour / number_of_neighbours;
        }
        mesh.vertices = updated_vertices;
    }
    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    std::cout << "Time elapsed: " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "ms" << std::endl;

    mesh.writeOBJ("goethe_out.obj");
}
#ifndef POLYGON_H
#define POLYGON_H

#define _CRT_SECURE_NO_WARNINGS 1

#include <vector>
#include <cstdio>
#include <string>
#include "vector.h"
#include <iostream>

// if the Polygon class name conflicts with a class in wingdi.h on Windows, use a namespace or change the name
class Polygon {  
public:
    double area()
    {
        if (vertices.size() < 3) return 0;
        double area_computed = 0;
        for (int i = 0; i < vertices.size(); i++)
        {
            area_computed += vertices[i][0] * vertices[(i + 1) % vertices.size()][1] - vertices[i][1] * vertices[(i + 1) % vertices.size()][0];
        }
        return area_computed / 2;
    }

    double integrate_squared_distance(const Vector& P)
    {
        double result = 0;
        if (vertices.size() < 3) return 0;
        for (int i = 1; i < vertices.size()-1; i++)
        {
            Vector triangle[3] = {vertices[0], vertices[i], vertices[i+1]};

            double local_result = 0;
            for (int k = 0; k < 3; k++)
            {
                for (int l = k; l < 3; l++)
                {
                    local_result += dot(triangle[k] - P, triangle[l] - P);
                }
            }

            double area_triangle = 0.5 * 
                abs(
                    (triangle[1][0] - triangle[0][0]) * (triangle[2][1] - triangle[0][1]) - 
                    (triangle[1][1] - triangle[0][1]) * (triangle[2][0] - triangle[0][0])
                );

            result += local_result / 6. * area_triangle;
        }
        return result;
    }

    Vector centroid()
    {
        Vector c = Vector(0, 0);
        int N = vertices.size();
        double A = area();
        for (int i = 0 ; i < N; i++)
        {
            c[0] += (vertices[i][0] + vertices[(i + 1) % N][0]) * (vertices[i][0] * vertices[(i + 1) % N][1] - vertices[(i + 1) % N][0] * vertices[i][1]);
            c[1] += (vertices[i][1] + vertices[(i + 1) % N][1]) * (vertices[i][0] * vertices[(i + 1) % N][1] - vertices[(i + 1) % N][0] * vertices[i][1]);
        }
        return c / (6 * A);
    }

    bool contains(const Vector& P)
    {
        int N = vertices.size();
        for (int i = 0; i < N; i++)
        {
            Vector v1 = vertices[(i + 1) % N];
            Vector v2 = P - vertices[i];
            if (v2[0] * (v1[1] - vertices[i][1]) > v2[1]*(vertices[i][0] - v1[0])) return false;
        }
        return true;
    }

	std::vector<Vector> vertices;
};	

// saves a static svg file. The polygon vertices are supposed to be in the range [0..1], and a canvas of size 1000x1000 is created
void save_svg(const std::vector<Polygon> &polygons, std::string filename, std::string fillcol = "none") {
    FILE* f = fopen(filename.c_str(), "w+"); 
    fprintf(f, "<svg xmlns = \"http://www.w3.org/2000/svg\" width = \"1000\" height = \"1000\">\n");
    for (int i=0; i<polygons.size(); i++) {
        fprintf(f, "<g>\n");
        fprintf(f, "<polygon points = \""); 
        for (int j = 0; j < polygons[i].vertices.size(); j++) {
            fprintf(f, "%3.3f, %3.3f ", (polygons[i].vertices[j][0] * 1000), (1000 - polygons[i].vertices[j][1] * 1000));
        }
        fprintf(f, "\"\nfill = \"%s\" stroke = \"black\"/>\n", fillcol.c_str());
        fprintf(f, "</g>\n");
    }
    fprintf(f, "</svg>\n");
    fclose(f);
}


// Adds one frame of an animated svg file. frameid is the frame number (between 0 and nbframes-1).
// polygons is a list of polygons, describing the current frame.
// The polygon vertices are supposed to be in the range [0..1], and a canvas of size 1000x1000 is created
void save_svg_animated(const std::vector<Polygon> &polygons, std::string filename, int frameid, int nbframes) {
    FILE* f;
    if (frameid == 0) {
        f = fopen(filename.c_str(), "w+");
        fprintf(f, "<svg xmlns = \"http://www.w3.org/2000/svg\" width = \"1000\" height = \"1000\">\n");
        fprintf(f, "<g>\n");
    } else {
        f = fopen(filename.c_str(), "a+");
    }
    fprintf(f, "<g>\n");
    for (int i = 0; i < polygons.size(); i++) {
        fprintf(f, "<polygon points = \""); 
        for (int j = 0; j < polygons[i].vertices.size(); j++) {
            fprintf(f, "%3.3f, %3.3f ", (polygons[i].vertices[j][0] * 1000), (1000-polygons[i].vertices[j][1] * 1000));
        }
        fprintf(f, "\"\nfill = \"none\" stroke = \"black\"/>\n");
    }
    fprintf(f, "<animate\n");
    fprintf(f, "	id = \"frame%u\"\n", frameid);
    fprintf(f, "	attributeName = \"display\"\n");
    fprintf(f, "	values = \"");
    for (int j = 0; j < nbframes; j++) {
        if (frameid == j) {
            fprintf(f, "inline");
        } else {
            fprintf(f, "none");
        }
        fprintf(f, ";");
    }
    fprintf(f, "none\"\n	keyTimes = \"");
    for (int j = 0; j < nbframes; j++) {
        fprintf(f, "%2.3f", j / (double)(nbframes));
        fprintf(f, ";");
    }
    fprintf(f, "1\"\n	dur = \"5s\"\n");
    fprintf(f, "	begin = \"0s\"\n");
    fprintf(f, "	repeatCount = \"indefinite\"/>\n");
    fprintf(f, "</g>\n");
    if (frameid == nbframes - 1) {
        fprintf(f, "</g>\n");
        fprintf(f, "</svg>\n");
    }
    fclose(f);
}

#endif
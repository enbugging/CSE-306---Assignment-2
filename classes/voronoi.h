#ifndef VORONOI_H
#define VORONOI_H

#define _CRT_SECURE_NO_WARNINGS 1

#include <string>
#include <vector>
#include "vector.h"
#include "polygon.h"

class Voronoi
{
public:
    Voronoi(const std::vector<Vector>& pts)
    {
        points = pts;
        for (Vector p : pts)
        {
            MIN_X = std::min(MIN_X, p[0]);
            MAX_X = std::max(MAX_X, p[0]);
            MIN_Y = std::min(MIN_Y, p[1]);
            MAX_Y = std::max(MAX_Y, p[1]);
        }
        MIN_X -= 1, MAX_X += 1, MIN_Y -= 1, MAX_Y += 1;
    }

    Polygon clip_by_bisector(const Polygon& poly, const Vector& P0, const Vector& Pi)
    {
        Polygon result;
        for (int i = 0; i < poly.vertices.size(); i++)
        {
            const Vector& A = (i == 0) ? poly.vertices[poly.vertices.size() - 1] : poly.vertices[i - 1];
            const Vector& B = poly.vertices[i];
            Vector M = (P0 + Pi)/2;
            double t = dot(M - A, Pi - P0)/dot(B - A, Pi - P0);
            Vector P = A + t*(B - A);

            if ((B - P0).norm() < (B - Pi).norm())
            {
                if ((A - P0).norm() > (A - Pi).norm()) result.vertices.push_back(P);
                result.vertices.push_back(B);
            }
            else if ((A - P0).norm() < (A - Pi).norm())
            {
                result.vertices.push_back(P);
            }
        }
        return result;
    }

    Polygon compute_voronoi_cell(int idx)
    {
        Polygon result;
        result.vertices.resize(4);
        result.vertices[0] = Vector(MIN_X, MIN_Y);
        result.vertices[1] = Vector(MIN_X, MAX_Y);
        result.vertices[2] = Vector(MAX_X, MAX_Y);
        result.vertices[3] = Vector(MAX_X, MIN_Y);

        for (int i = 0; i < points.size(); i++)
        {
            if (idx == i) continue;
            result = clip_by_bisector(result, points[idx], points[i]);
        }

        return result;
    }

    void compute()
    {
        voronoi.resize(points.size());
        for (int i = 0; i < points.size(); i++)
        {
            voronoi[i] = compute_voronoi_cell(i);
        }
    }

    void save(std::string filename)
    {
        save_svg(voronoi, filename, "blue");
    }

    std::vector<Vector> points;
    std::vector<Polygon> voronoi;
    double MIN_X, MAX_X, MIN_Y, MAX_Y;
};

#endif
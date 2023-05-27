#ifndef VORONOI_H
#define VORONOI_H

#define _CRT_SECURE_NO_WARNINGS 1

#include <string>
#include <vector>
#include <algorithm>
#include <limits>
#include <utility>
#include "vector.h"
#include "polygon.h"
#include "kd_tree.h"

class Voronoi
{
public:
    Voronoi(int N_disc = 50): boundary_computed(false) {
        disc.vertices.resize(N_disc);
        for (int i = 0; i < N_disc; i++)
        {
            double theta = 2*M_PI*i/N_disc;
            disc.vertices[i] = Vector(cos(theta), sin(theta));
        }
    }

    Voronoi(const std::vector<Vector>& pts, const std::vector<double>& weights, int N_disc = 50)
    {
        points = pts;
        this->weights = weights;
        this->boundary_computed = true;

        disc.vertices.resize(N_disc);
        for (int i = 0; i < N_disc; i++)
        {
            double theta = 2*M_PI*i/N_disc;
            disc.vertices[i] = Vector(cos(theta), sin(theta));
        }

        MIN_X = std::numeric_limits<double>::max();
        MAX_X = std::numeric_limits<double>::min();
        MIN_Y = std::numeric_limits<double>::max();
        MAX_Y = std::numeric_limits<double>::min();
        for (Vector p : pts)
        {
            MIN_X = std::min(MIN_X, p[0]);
            MAX_X = std::max(MAX_X, p[0]);
            MIN_Y = std::min(MIN_Y, p[1]);
            MAX_Y = std::max(MAX_Y, p[1]);
        }
        MIN_X -= 1e-4, MAX_X += 1e-4, MIN_Y -= 1e-4, MAX_Y += 1e-4;

        max_weight = (*std::max_element(weights.begin(), weights.end()));
        std::vector<std::pair<Vector, int> > aux_points;
        aux_points.resize(points.size());
        for(int i = 0; i < points.size(); i++)
        {
            aux_points[i] = std::make_pair(points[i], i);
            aux_points[i].first[2] = sqrt(max_weight - weights[i]);
        }
        kdtree = KDTree(aux_points);
    }

    Polygon clip_by_bisector(
        const Polygon& poly, 
        int idx_0, 
        int idx_i, 
        const Vector& P0, 
        const Vector& Pi, 
        double& maximum_distance)
    {
        Polygon result;
        Vector M = (P0 + Pi)/2;
        Vector Mp = M + (weights[idx_0] - weights[idx_i])/(2.*(P0 - Pi).norm2())*(Pi - P0);
        
        maximum_distance = 0;
        for (int i = 0; i < poly.vertices.size(); i++)
        {
            const Vector& A = (i == 0) ? poly.vertices[poly.vertices.size() - 1] : poly.vertices[i - 1];
            const Vector& B = poly.vertices[i];
            
            double t = dot(Mp - A, Pi - P0)/dot(B - A, Pi - P0);
            Vector P = A + t*(B - A);

            if ((B - P0).norm2() - weights[idx_0] < (B - Pi).norm2() - weights[idx_i])
            {
                if ((A - P0).norm2() - weights[idx_0] > (A - Pi).norm2() - weights[idx_i])
                {
                    result.vertices.push_back(P);
                    maximum_distance = std::max(maximum_distance, (P - P0).norm2());
                }
                result.vertices.push_back(B);
                maximum_distance = std::max(maximum_distance, (B - P0).norm2());
            }
            else if ((A - P0).norm2() - weights[idx_0] < (A - Pi).norm2() - weights[idx_i])
            {
                result.vertices.push_back(P);
                maximum_distance = std::max(maximum_distance, (P - P0).norm2());
            }
        }
        //for (Vector& p : result.vertices)
        //{
        //    std::cout << p[0] << " " << p[1] << " - " << P0[0] << " " << P0[1] << " " << (p - P0).norm2() << std::endl;
        //}
        maximum_distance += max_weight - weights[idx_0];
        return result;
    }

    Polygon clip_by_edge(const Polygon& poly, const Vector& u, const Vector& v) const
    {
        Polygon result;
        Vector N(v[1] - u[1], u[0] - v[0]);
        for (int i = 0; i < poly.vertices.size(); i++)
        {
            const Vector& A = (i == 0) ? poly.vertices[poly.vertices.size() - 1] : poly.vertices[i - 1];
            const Vector& B = poly.vertices[i];
            
            double t = dot(u - A, N)/dot(B - A, N);
            Vector P = A + t*(B - A);

            if (dot(B - u, N) < 0) // B is inside
            {
                if (dot(A - u, N) > 0) result.vertices.push_back(P); // A is outside
                result.vertices.push_back(B);
            }
            else if (dot(A - u, N) < 0) // A is inside
            {
                result.vertices.push_back(P);
            }
        }
        return result;
    }

    Polygon intersect_with_disc(const Polygon& polygon, const Vector center, double radius) const
    {
        Polygon result(polygon);
        for (int i = 0; i < disc.vertices.size(); i++)
        {
            const Vector& u = disc.vertices[i];
            const Vector& v = i == disc.vertices.size() - 1 ? disc.vertices[0] : disc.vertices[i + 1];
            result = clip_by_edge(result, center + radius*u, center + radius*v);
        }
        return result;
    }

    Polygon compute_voronoi_cell(int idx)
    {
        Polygon result;
        double maximum_distance = std::numeric_limits<double>::max();
        bool contributing = true;
        result.vertices.resize(4);
        result.vertices[0] = Vector(MIN_X, MIN_Y);
        result.vertices[1] = Vector(MAX_X, MIN_Y);
        result.vertices[2] = Vector(MAX_X, MAX_Y);
        result.vertices[3] = Vector(MIN_X, MAX_Y);
        /* Version without kd-tree
        for (int i = 0; i < points.size(); i++)
        {
            if (idx == i) continue;
            result = clip_by_bisector(result, idx, i, points[idx], points[i], maximum_distance);
        }
        //*/
        //* Version with kd-tree
        const int k = 50;
        std::vector<int> points_in_range;
        //::cout << "building for " << idx << "\n";
        double distance = 0;
        for (int prev_k = 0, curr_k = k; prev_k < points.size() - 1; curr_k = prev_k + k)
        {
            //std::cout << "retrieving for " << curr_k << "\n";
            contributing = true;
            kdtree.findKNearestNeighbors(points_in_range, points[idx], curr_k);

            for (; prev_k < points_in_range.size(); prev_k++)
            {
                if (idx == points_in_range[prev_k]) continue;

                double d = sqrt(max_weight - weights[points_in_range[prev_k]]) - sqrt(max_weight - weights[idx]);
                d *= d;
                //std::cout << d << " " << (points[idx] - points[points_in_range[prev_k]]).norm2() << "\n";
                //std::cout << (d + (points[idx] - points[points_in_range[prev_k]]).norm2()) << " " << 2*maximum_distance << "\n";
                if (d + (points[idx] - points[points_in_range[prev_k]]).norm2() > 2*maximum_distance)
                {
                    contributing = false;
                    break;
                }

                result = clip_by_bisector(result, idx, points_in_range[prev_k], points[idx], points[points_in_range[prev_k]], maximum_distance);
            }
            if (not contributing) break;
        }
        //*/
        //result = intersect_with_disc(result, points[idx], sqrt(weights[idx] - weights[weights.size() - 1]));
        return result;
    }

    void compute()
    {
        if (not boundary_computed)
        {
            this -> boundary_computed = true;
            MIN_X = std::numeric_limits<double>::max();
            MAX_X = std::numeric_limits<double>::min();
            MIN_Y = std::numeric_limits<double>::max();
            MAX_Y = std::numeric_limits<double>::min();
            for (Vector p : points)
            {
                MIN_X = std::min(MIN_X, p[0]);
                MAX_X = std::max(MAX_X, p[0]);
                MIN_Y = std::min(MIN_Y, p[1]);
                MAX_Y = std::max(MAX_Y, p[1]);
            }
            MIN_X -= 1e-4, MAX_X += 1e-4, MIN_Y -= 1e-4, MAX_Y += 1e-4;

            max_weight = (*std::max_element(weights.begin(), weights.end()));
            std::vector<std::pair<Vector, int> > aux_points;
            aux_points.resize(points.size());
            for(int i = 0; i < points.size(); i++)
            {
                aux_points[i] = std::make_pair(points[i], i);
                aux_points[i].first[2] = sqrt(max_weight - weights[i]);
            }
            kdtree = KDTree(aux_points);
        }
        voronoi.resize(points.size());
        #pragma omp parallel for schedule(dynamic)
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
    std::vector<double> weights;
    double MIN_X, MAX_X, MIN_Y, MAX_Y;
    double max_weight;
    bool boundary_computed;
    Polygon disc;
    KDTree kdtree;
};

#endif
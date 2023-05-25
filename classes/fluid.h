#ifndef FLUID_H
#define FLUID_H

#define STB_IMAGE_WRITE_IMPLEMENTATION

#include <vector>
#include <sstream>
#include <string>

#include "vector.h"
#include "polygon.h"
#include "optimal_transport.h"
#include "stb_image_write.h"

int sgn(double x) {
    if (x > 0) return 1;
    if (x < 0) return -1;
    return 0;
}

void save_frame(const std::vector<Polygon> &cells, std::string filename, int frameid = 0, int N = 60) {
    int W = 500, H = 500;
    std::vector<unsigned char> image(W*H * 3, 255);
#pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < cells.size(); i++) {

        double bminx = 1E9, bminy = 1E9, bmaxx = -1E9, bmaxy = -1E9;
        for (int j = 0; j < cells[i].vertices.size(); j++) {
            bminx = std::min(bminx, cells[i].vertices[j][0]);
            bminy = std::min(bminy, cells[i].vertices[j][1]);
            bmaxx = std::max(bmaxx, cells[i].vertices[j][0]);
            bmaxy = std::max(bmaxy, cells[i].vertices[j][1]);
        }
        bminx = std::min(W-1., std::max(0., W * bminx));
        bminy = std::min(H-1., std::max(0., H * bminy));
        bmaxx = std::max(W-1., std::max(0., W * bmaxx));
        bmaxy = std::max(H-1., std::max(0., H * bmaxy));

        for (int y = bminy; y < bmaxy; y++) {
            for (int x = bminx; x < bmaxx; x++) {
                int prevSign = 0;
                bool isInside = true;
                double mindistEdge = 1E9;
                for (int j = 0; j < cells[i].vertices.size(); j++) {
                    double x0 = cells[i].vertices[j][0] * W;
                    double y0 = cells[i].vertices[j][1] * H;
                    double x1 = cells[i].vertices[(j + 1) % cells[i].vertices.size()][0] * W;
                    double y1 = cells[i].vertices[(j + 1) % cells[i].vertices.size()][1] * H;
                    double det = (x - x0)*(y1-y0) - (y - y0)*(x1-x0);
                    int sign = sgn(det);
                    if (prevSign == 0) prevSign = sign; else
                        if (sign == 0) sign = prevSign; else
                        if (sign != prevSign) {
                            isInside = false;
                            break;
                        }
                    prevSign = sign;
                    double edgeLen = sqrt((x1 - x0)*(x1 - x0) + (y1 - y0)*(y1 - y0));
                    double distEdge = std::abs(det)/ edgeLen;
                    double dotp = (x - x0)*(x1 - x0) + (y - y0)*(y1 - y0);
                    if (dotp<0 || dotp>edgeLen*edgeLen) distEdge = 1E9;
                    mindistEdge = std::min(mindistEdge, distEdge);
                }
                if (isInside) {
                    if (i < N) {   // the N first particles may represent fluid, displayed in blue
                    	image[((H - y - 1)*W + x) * 3] = 0;
                    	image[((H - y - 1)*W + x) * 3 + 1] = 0;
                    	image[((H - y - 1)*W + x) * 3 + 2] = 255;
                    }
                    //if (mindistEdge <= 2) {
                    //    image[((H - y - 1)*W + x) * 3] = 0;
                    //    image[((H - y - 1)*W + x) * 3 + 1] = 0;
                    //    image[((H - y - 1)*W + x) * 3 + 2] = 0;
                    //}

                }
                
            }
        }
    }
    std::ostringstream os;
    os << filename << frameid << ".png";
    stbi_write_png(os.str().c_str(), W, H, 3, &image[0], 0);
}

class Fluid
{
public:
    Fluid() {}
    Fluid(int N, double MIN_X = 0, double MAX_X = 1, double MIN_Y = 0, double MAX_Y = 1): MIN_X(MIN_X), MAX_X(MAX_X), MIN_Y(MIN_Y), MAX_Y(MAX_Y)
    {
        particles.resize(N);
        velocities.resize(N);
        for (int i = 0; i < N; i++)
        {
            particles[i] = Vector(rand()/(double) RAND_MAX * (MAX_X - MIN_X) + MIN_X, rand()/(double) RAND_MAX * (MAX_Y - MIN_Y) + MIN_Y);
            velocities[i] = Vector(0, 0);
        }
    }

    void stepfluid()
    {
        otsolver.pts = particles;
        otsolver.lambdas = std::vector<double>(particles.size(), 1./particles.size() * VOLUME_FLUID);
        otsolver.solve();

        const double 
            mass_particles = 200, 
            epsilon = 0.004, 
            epsilon2 = epsilon*epsilon, 
            dt = 0.002;
        Vector next_particle;
        for (int i = 0; i < particles.size(); i++)
        {
            Vector gravity = Vector(0, -9.80665) * mass_particles;
            Vector centroid = otsolver.solution.voronoi[i].centroid();
            Vector otForces = 1./epsilon2 * (centroid - particles[i]) * otsolver.solution.weights[i] * mass_particles;
            Vector forces = gravity + otForces;
            velocities[i] += dt / mass_particles * forces;
            next_particle = particles[i] + dt * velocities[i];
            // bounce particle back into domain if it reaches the boundary
            if (next_particle.data[0] < MIN_X ||
                next_particle.data[0] > MAX_X ||
                next_particle.data[1] < MIN_Y ||
                next_particle.data[1] > MAX_Y) velocities[i] = -velocities[i];
            //std::cout << "particle " << i << " " << next_particle.data[0] << " " << next_particle.data[1] << std::endl;
            particles[i] += dt * velocities[i];
            //std::cout << "particle " << i << " " << particles[i].data[0] << " " << particles[i].data[1] << std::endl;
        }
    }

    void runfluid(std::string filename = "fluid/animation")
    {
        for (int i = 0; i < 200; i++)
        {
            //std::cout << "Iteration: " << i << std::endl;
            stepfluid();
            save_frame(otsolver.solution.voronoi, filename, i);
        }
    }

    OptimalTransport otsolver;
    std::vector<Vector> particles;
    std::vector<Vector> velocities;
    double MIN_X, MAX_X, MIN_Y, MAX_Y;
};

#endif
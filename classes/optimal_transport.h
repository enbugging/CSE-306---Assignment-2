#ifndef OPTIMAL_TRANSPORT_H
#define OPTIMAL_TRANSPORT_H

#include <vector>
#include "vector.h"
#include "voronoi.h"
#include "lbfgs.h"
#include <iostream>

#define VOLUME_AIR 0.7
#define VOLUME_FLUID 0.3

class OptimalTransport {
public:
    OptimalTransport() {}

    OptimalTransport(const std::vector<Vector>& pts, const std::vector<double>& lambdas) {
        this->pts = pts;
        this->lambdas = lambdas;
    }

    static lbfgsfloatval_t _evaluate(
        void *instance,
        const lbfgsfloatval_t *x,
        lbfgsfloatval_t *g,
        const int n,
        const lbfgsfloatval_t step
        )
    {
        return reinterpret_cast<OptimalTransport*>(instance)->evaluate(x, g, n, step);
    }

    lbfgsfloatval_t evaluate(
        const lbfgsfloatval_t *x,
        lbfgsfloatval_t *g,
        const int n,
        const lbfgsfloatval_t step
        )
    {
        // Compute the function value
        lbfgsfloatval_t fx = 0.0;

        for (int i = 0; i < n; i++ )
        {
            solution.weights[i] = x[i];
        }
        solution.compute();

        double s1 = 0.0, s2 = 0.0, s3 = 0.0, estimated_volume_fluid = 0.0;
        for (int i = 0; i < n-1; i++)
        {
            //std::cout << "Polygon size: " << solution.voronoi[i].vertices.size() << '\n';
            double A = solution.voronoi[i].area();
            g[i] = -(lambdas[i] - A);
            s1 += solution.voronoi[i].integrate_squared_distance(solution.points[i]);
            s2 -= A * x[i];
            //std::cout << A << '\n';
            s3 += lambdas[i] * x[i];
            estimated_volume_fluid += A;
        }
        fx = s1 + s2 + s3;
        // account for volume of air, i.e. w_air (desired_volume - esimated_volume)
        fx += x[n-1] * (VOLUME_AIR - (1-estimated_volume_fluid));
        g[n-1] = -(VOLUME_AIR - (1-estimated_volume_fluid));
        return -fx;
    }

    static int _progress(
        void *instance,
        const lbfgsfloatval_t *x,
        const lbfgsfloatval_t *g,
        const lbfgsfloatval_t fx,
        const lbfgsfloatval_t xnorm,
        const lbfgsfloatval_t gnorm,
        const lbfgsfloatval_t step,
        int n,
        int k,
        int ls
        )
    {
        return reinterpret_cast<OptimalTransport*>(instance)->progress(x, g, fx, xnorm, gnorm, step, n, k, ls);
    }

    int progress(
        const lbfgsfloatval_t *x,
        const lbfgsfloatval_t *g,
        const lbfgsfloatval_t fx,
        const lbfgsfloatval_t xnorm,
        const lbfgsfloatval_t gnorm,
        const lbfgsfloatval_t step,
        int n,
        int k,
        int ls
        )
    {
        //printf("Iteration %d:\n", k);
        //printf("  fx = %f, x[0] = %f, x[1] = %f\n", fx, x[0], x[1]);
        //printf("  xnorm = %f, gnorm = %f, step = %f\n", xnorm, gnorm, step);
        //printf("\n");
        return 0;
    }

    void solve() {
        solution.points = pts;
        solution.weights.resize(pts.size()+1);
        std::fill(solution.weights.begin(), solution.weights.end(), 1.0);
        solution.weights[pts.size()] = 0.999;

        solution.compute();
        double fx = 0.0;

        // L-BFGS
        size_t ret = lbfgs((int) pts.size() + 1, &solution.weights[0], &fx, _evaluate, _progress, this, NULL);
    }

    std::vector<Vector> pts;
    std::vector<double> lambdas;
    Voronoi solution;
};

#endif
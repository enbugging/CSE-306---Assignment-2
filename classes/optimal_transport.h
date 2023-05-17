#ifndef OPTIMAL_TRANSPORT_H
#define OPTIMAL_TRANSPORT_H

#include <vector>
#include "vector.h"
#include "voronoi.h"
#include "lbfgs.h"
#include <iostream>

class OptimalTransport {
public:
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

        double s1 = 0.0, s2 = 0.0, s3 = 0.0;
        for (int i = 0; i < n; i++)
        {
            //std::cout << "Polygon size: " << solution.voronoi[i].vertices.size() << '\n';
            double A = solution.voronoi[i].area();
            g[i] = -(lambdas[i] - A);
            s1 += solution.voronoi[i].integrate_squared_distance(solution.points[i]);
            s2 -= A * x[i];
            //std::cout << A << '\n';
            s3 += lambdas[i] * x[i];
        }
        fx = s1 + s2 + s3;
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
        solution.weights.resize(pts.size());
        std::fill(solution.weights.begin(), solution.weights.end(), 1.0);
        std::cout << solution.points.size() << '\n';
        solution.compute();
        std::cout << solution.voronoi[0].vertices.size() << '\n';
        double fx = 0.0;

        // L-BFGS
        size_t ret = lbfgs((int) pts.size(), &solution.weights[0], &fx, _evaluate, _progress, this, NULL);

        solution.compute();
    }

    std::vector<Vector> pts;
    std::vector<double> lambdas;
    Voronoi solution;
};

#endif
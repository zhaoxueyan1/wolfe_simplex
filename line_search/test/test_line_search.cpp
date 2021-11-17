//
// Created by zhaox on 2021/11/17.
//
#include "Eigen/Dense"
#include <iostream>

namespace GoldenSearch {
    using namespace Eigen;

    VectorXd f(VectorXd x) {
        VectorXd res(1);
        res << -std::min(std::min(x[0] / 2, 2 - pow(x[0] - 3, 2)), 2 - x[0] / 2);
        return res;
    }

    VectorXd nebula_f(VectorXd x) {
        double delta = 0.001;
        VectorXd t(1);
        t << delta;
        return (f(x + t) - f(x)) / delta;
    }

    double get_fix_step_alpha(VectorXd x, VectorXd d) {
        double delta_s = 0.5;
        int i = 1;
        for (; i < 10; i++) {
            if (f(x + delta_s * i * d)[0] < f(x + delta_s * (i + 1) * d)[0]) {
                break;
            }
        }
        return i * delta_s;
    }

    VectorXd line_search_solve(VectorXd x, int itr) {
        for (int i = 0; i < itr; i++) {
            VectorXd d = -nebula_f(x);
            double alpha = get_fix_step_alpha(x, d);
            std::cout << "itr: " << i << std::endl;
            std::cout << "d:" << d << std::endl;
            std::cout << "alpha:" << alpha << std::endl;
            std::cout << "X:" << x << std::endl;
            printf("F:%f\n", f(x)[0]);
            x = x + alpha * d;
            x[0] = std::max(x[0], 0.0);
            x[0] = std::min(x[0], 8.0);
        }
        return x;
    }

    double golden_search_solve(double l, double r) {

//        printf("%f %f\n", l, r);
        VectorXd x(1);
        double b1 = 0.618 * l + (1 - 0.618) * r;
        double c1 = (1 - 0.618) * l + 0.618 * r;
        x << b1;
        double f1 = f(x)[0];
        x << c1;
        double f2 = f(x)[0];
        while (r - l > 1e-3) {
            printf("l:%f r:%f f1: %f f2: %f\n", l, r, f1, f2);
            if (f1 <= f2) {
                r = c1;
                f2 = f1;
                c1 = b1;
                b1 = 0.618 * l + (1 - 0.618) * r;
                x << b1;
                f1 = f(x)[0];
                printf("f1 %f\n", f1);
            } else {
                l = b1;
                f1 = f2;
                b1 = c1;
                c1 = (1 - 0.618) * l + 0.618 * r;
                x << c1;
                f2 = f(x)[0];
                printf("f2 %f\n", f2);
            }
        }
        return r;
    }

    int fib[100];

    void gen_fib(int n) {
        fib[0] = 0;
        fib[1] = 1;
        for (int i = 2; i <= n; i++) {
            fib[i] = fib[i - 1] + fib[i - 2];
        }
    }

    double fibonacci_search_solve(double l, double r) {
        int n = 8;
        gen_fib(n);
//        printf("%f %f\n", l, r);
        VectorXd x(1);
        int itr = 0;
        double lambda = 1. * fib[n - itr - 1] / fib[n - itr];
        double b1 = lambda * l + (1 - lambda) * r;
        double c1 = (1 - lambda) * l + lambda * r;
        x << b1;
        double f1 = f(x)[0];
        x << c1;
        double f2 = f(x)[0];
        while (r - l > 0.1) {
            printf("itr: %d\nl:%f b:%f c:%f r:%f f1: %f f2: %f\n", itr, l, b1, c1, r, f1, f2);
            itr++;
            lambda = 1. * fib[n - itr - 1] / fib[n - itr];
            if (f1 <= f2) {
                r = c1;
                f2 = f1;
                c1 = b1;
                b1 = lambda * l + (1 - lambda) * r;
                x << b1;
                f1 = f(x)[0];
                printf("f1 %f\n", f1);
            } else {
                l = b1;
                f1 = f2;
                b1 = c1;
                c1 = (1 - lambda) * l + lambda * r;
                x << c1;
                f2 = f(x)[0];
                printf("f2 %f\n", f2);
            }

            if (c1 - b1 < 1e-7) {
                break;
            }
        }
        return r;
    }
}

namespace Zigzag {
    using namespace Eigen;

    VectorXd f(VectorXd x) {
        VectorXd res(1);
        res << pow(x[0] - 2, 2) + 4 * pow(x[1] - 3, 2);
        return res;
    }

    VectorXd nebula_f(VectorXd x) {
        VectorXd res(2);
        res << 2 * (x[0] - 2), 8 * (x[1] - 3);
        return res;
    }

    VectorXd line_search_solve(VectorXd x, int itr) {
        for (int i = 0; i < itr; i++) {
            VectorXd n_f = nebula_f(x);
            double norm = sqrt((n_f.transpose() * n_f)[0]);
            VectorXd d = -n_f / norm;
            double alpha = 1;
            printf("|%d |", i);
            printf("(%f , %f) |", x[0], x[1]);
            printf("%f |", f(x)[0]);
            printf("(%f,%f) |", n_f[0], n_f[1]);
            printf("%f |\n", norm);

            x = x + alpha * d;
        }
        return x;
    }

}

# if 0
Eigen::VectorXd x(2);
    x << 0, 0;
    Zigzag::line_search_solve(x, 10);
//    GoldenSearch::fibonacci_search_solve(0, 8);
#endif
//
// Created by zhaox on 2021/11/17.
//
#include "Eigen/Dense"
#include "../../qp/interior_point_qp/QuadraticProgramming.h"
namespace SQP {
    using namespace Eigen;

    VectorXd f(VectorXd x) {
        VectorXd res(1);
        res << pow(x[0] - 9.0 / 4.0, 2) + pow(x[1] - 2, 2);
        return res;
    }

    VectorXd nebula_f(VectorXd x) {
        VectorXd res(2);
        res << 2 * (x[0] - 9.0 / 4.0), 2 * (x[1] - 2);
        return res;
    }

    VectorXd g(VectorXd x) {
        VectorXd res(4);
        res << (x[0] + x[1] - 6),
                (x[0] * x[0] - x[1]),
                (-x[0]),
                (-x[1]);
        return res;
    }

    MatrixXd nebula_g(VectorXd x) {
        MatrixXd res(4, 2);
        res << 1.0, 1.0,
                2 * x[0], -1.0,
                -1.0, 0,
                0, -1.0;
        return res;
    }

    VectorXd f_p(VectorXd x) {
        VectorXd res(1);
        VectorXd t_g = g(x);
        double Max = std::max(std::max(t_g[0], t_g[1]), std::max(t_g[2], t_g[3]));
        res << f(x)[0] + (Max > 1e-7 ? 150000.0 : 0.0);
        return res;
    }


    VectorXd L(VectorXd x, VectorXd lambda) {
        VectorXd res(1);
        res << f(x)[0] + lambda[1] * (x[0] * x[0] - x[1]) + lambda[0] * (x[0] + x[1] - 6) + lambda[2] * (6 - x[0]) +
               lambda[3] * (6 - x[1]);
        return res;
    }

    VectorXd nebula_x_L(VectorXd x, VectorXd lambda) {
        VectorXd res(2);
        res << 2 * (x[0] - 9.0 / 4.0) + lambda[0] + 2 * lambda[1] * x[0] - lambda[2], 2 * (x[1] - 2) + lambda[0] -
                                                                                      lambda[1] - lambda[3];
        return res;
    }

    MatrixXd hessian_x_L(VectorXd x, VectorXd lambda) {
        MatrixXd res(2, 2);
        res << 2 + 2 * lambda[1], 0,
                0, 2;
        return res;
    }

    double get_fix_step_alpha(VectorXd x, VectorXd d) {
        double delta_s = 0.1;
        int i = 1;
        for (; i < 10; i++) {
            VectorXd x0 = x + delta_s * i * d;
            if (f_p(x0)[0] < f_p(x + delta_s * (i + 1) * d)[0]) {
                break;
            }
        }
        return i * delta_s;
    }


    VectorXd solve(int itr) {
        VectorXd x(2);
        x << 0, 0;
        VectorXd lambda(4);
        lambda << 0, 0, 0, 0;
        for (int k = 0; k < itr; k++) {
            MatrixXd H_k = hessian_x_L(x, lambda);
            VectorXd nebula_k = nebula_f(x);
            MatrixXd n_g_k = nebula_g(x);
            VectorXd g_k = g(x);
            MatrixXd A(0, 0);
            VectorXd b(0);

            QuadraticProgramming qp(H_k, nebula_k, n_g_k, -g_k, A, b);
            VectorXd init_x(2);
            init_x << 0, 0;
            VectorXd d = qp.solve(init_x);
//            std::cout<<"Constraint\n"<<n_g_k * d+ g_k<<std::endl;
            VectorXd x_optimal = x + d;
            VectorXd temp = -H_k * x_optimal - nebula_k;

            VectorXd lambda_optimal(4);
//            std::cout<<n_g_k << std::endl;
            lambda_optimal = n_g_k.transpose().bdcSvd(ComputeThinU | ComputeThinV).solve(temp);
            double alpha = get_fix_step_alpha(x, d);
            printf("|%d |", k);
            printf("(%f , %f) |", x[0], x[1]);
            printf("%f |", f(x)[0]);
            printf("(%f,%f) |", d[0], d[1]);
            printf("%f |\n", alpha);

            x = x + alpha * d;
            lambda = lambda + alpha * (lambda_optimal - lambda);

        }
        return x;
    }
}

namespace convex_interior_point {
    using namespace Eigen;

    VectorXd f(VectorXd x) {
        VectorXd res(1);
        res << -x[0] * x[1];
        return res;
    }

    VectorXd nebula_f(VectorXd x) {
        VectorXd res(2);
        res << 2 * (x[0] - 9.0 / 4.0), 2 * (x[1] - 2);
        return res;
    }

    VectorXd g(VectorXd x) {
        VectorXd res(1);
        res << (-1 + x[0] * x[0] + x[1] * x[1]);
        return res;
    }

}

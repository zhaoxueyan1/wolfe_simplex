
#include <iostream>
#include "wolfe_simplex/wolfe_simplex.h"
#include "active_set_dual/QuadProg.h"
#include "interior_point_qp/QuadraticProgramming.h"
//#include "opt_unconstr/newton.h"

void test_wolfe_simplex_method() {
    SHG::Vecdouble p((std::vector<double>) {-2, -5, 0, 0});
    SHG::Vecdouble C((std::vector<double>) {1, 0, 0, 0,
                                            1, 0, 0,
                                            0, 0,
                                            0});
    SHG::Matdouble A(2, 4, SHG::Vecdouble((std::vector<double>) {1, -2, -1, 0,
                                                                 -1, -2, 0, -1,
    }));
    SHG::Vecdouble b((std::vector<double>) {-2, -2});
    SHG::Vecdouble x(4);
    double f = 0;
    int ok = SHG::wolfe(p, C, A, b, x, f);
    if (!ok) {
        printf("OK! The Wolfe Simplex method result is %f\n", f);
        for (int i = 0; i < x.size(); i++) {
            printf("X%d : %f ", i, x[i]);
        }
    } else {
        printf("Error %d\n", ok);
    }
    putchar(10);
}

void test_active_set_dual_method() {
    SHG::Vecdouble p((std::vector<double>) {-2, -5});
    SHG::Matdouble C(2, 2, SHG::Vecdouble((std::vector<double>{
            2, 0,
            0, 2})));
    SHG::Matdouble CI(4, 2, SHG::Vecdouble((std::vector<double>) {1, -2,
                                                                  -1, -2,
                                                                  1, 0,
                                                                  0, 1
    }));
    CI = transpose(CI);
    SHG::Vecdouble ci0((std::vector<double>) {2, 2, 0, 0});
    SHG::Matdouble CE(2, 0, SHG::Vecdouble((std::vector<double>) {0, 0}));
    SHG::Vecdouble b((std::vector<double>(0)));
    SHG::Vecdouble x(4);

    double f = SHG::solve_quadprog(C, p, CE, b, CI, ci0, x);
    printf("active-set dual method OK! The result is %f\n", f);
    for (int i = 0; i < x.size(); i++) {
        printf("X%d : %f ", i, x[i]);
    }
    putchar(10);
}

void test_primary_dual_interior_point_method() {
    using namespace Eigen;
    VectorXd p(2);
    p << -2, -5;
    MatrixXd C(2, 2);
    C << 2, 0,
            0, 2;
    MatrixXd G(4, 2);
    G << -1, 2,
            1, 2,
            -1, 0,
            0, -1;
    VectorXd h(4);
    h << 2, 2, 0, 0;

    MatrixXd A(0, 0);
    VectorXd b(0);

    VectorXd x = VectorXd::Zero(2);
    VectorXd f(2);
    std::cout << C << std::endl;
    std::cout << G << std::endl;
    QuadraticProgramming qp(C, p, G, h, A, b);
    f = qp.solve(x);
    std::cout << x << std::endl;
    VectorXd res_val(1);
    res_val = 0.5 * f.transpose() * C * f + p.transpose() * f;
    printf("OK! The result is %f\n", res_val[0]);
    for (int i = 0; i < x.size(); i++) {
        printf("X%d : %f ", i, f[i]);
    }
}

namespace Newton {
    using namespace Eigen;

    VectorXd f(VectorXd x) {
        VectorXd res(1);
        res << pow(x[0] - 3, 4) + pow(x[0] - 3 * x[1], 2);
        return res;
    }

    VectorXd nebula_f(VectorXd x) {
        VectorXd res(2);
        res << 4 * pow(x[0] - 3, 3) + 2 * (x[0] - 3 * x[1]), 6 * (3 * x[1] - x[0]);
        return res;
    }

    MatrixXd hessian_f(VectorXd x) {
        MatrixXd res(2, 2);
        double x1x1 = 12 * pow(x[0] - 3, 2) + 2;
        double x1x2 = -6;
        double x2x2 = 18;

        res << x1x1, x1x2,
                x1x2, x2x2;

        return res;
    }

    VectorXd newton_solve(VectorXd init_x, int itr) {
        Eigen::VectorXd res(init_x);
        for (int i = 0; i < itr; ++i) {
            VectorXd f_v = f(res);
            if (f_v[0] < 1e-7) {
                break;
            }
            VectorXd nf_v = nebula_f(res);
            MatrixXd H = hessian_f(res);
            std::cout << "itr:\n " << i << std::endl;
            std::cout << "f: \n" << f_v << std::endl;
            std::cout << "nebula_f: \n" << nf_v << std::endl;
            std::cout << "Hessian_f: \n" << H << std::endl;
            std::cout << "x: \n" << res << std::endl;
            res = res - H.inverse() * nf_v;
        }
        return res;
    }

    VectorXd Levenberg_Marquardt_solve(VectorXd init_x, int itr) {
        Eigen::VectorXd res(init_x);
        for (int i = 0; i < itr; ++i) {
            VectorXd f_v = f(res);
            if (f_v[0] < 1e-7) {
                break;
            }
            VectorXd nf_v = nebula_f(res);
            MatrixXd H = hessian_f(res);
            std::cout << "itr:\n " << i << std::endl;
            std::cout << "f: \n" << f_v << std::endl;
            std::cout << "nebula_f: \n" << nf_v << std::endl;
            std::cout << "Hessian_f: \n" << H << std::endl;
            std::cout << "x: \n" << res << std::endl;
            MatrixXd I(2, 2);
            I << 1, 0,
                    0, 1;
            double lambda = 1;
            res = res - (H + lambda * I).inverse() * nf_v;
        }
        return res;
    }

    VectorXd BFGS_solve(VectorXd init_x, MatrixXd H, int itr) {

        Eigen::VectorXd x(init_x);
        VectorXd s(2);
        VectorXd q(2);
        VectorXd x_new;
        for (int i = 0; i < itr; ++i) {
            VectorXd nf_v = nebula_f(x);
            VectorXd f_v = f(x);
            x_new = x - H.inverse() * nf_v;
            s = x_new - x;
            q = nebula_f(x_new) - nf_v;

            std::cout << "itr:\n " << i << std::endl;
            std::cout << "f: \n" << f_v << std::endl;
            std::cout << "nebula_f: \n" << nf_v << std::endl;
            std::cout << "Hessian_f: \n" << H << std::endl;
            std::cout << "x: \n" << x << std::endl;
            H = H + q * q.transpose() / (q.transpose() * s) -
                H * s * s.transpose() * H.transpose() / (s.transpose() * H * s);
            x = x_new;

        }
        return x;
    }

    VectorXd DFP_solve(VectorXd init_x, MatrixXd H, int itr) {

        Eigen::VectorXd x(init_x);
        VectorXd s(2);
        VectorXd q(2);
        VectorXd x_new;
        for (int i = 0; i < itr; ++i) {
            VectorXd nf_v = nebula_f(x);
            VectorXd f_v = f(x);
            x_new = x - H * nf_v;
            s = x_new - x;
            q = nebula_f(x_new) - nf_v;

            std::cout << "itr:\n " << i << std::endl;
            std::cout << "f: \n" << f_v << std::endl;
            std::cout << "nebula_f: \n" << nf_v << std::endl;
            std::cout << "Hessian_f: \n" << H.inverse() << std::endl;
            std::cout << "x: \n" << x << std::endl;
            H = H + s * s.transpose() / (q.transpose() * s) -
                H * q * q.transpose() * H.transpose() / (q.transpose() * H * q);
            x = x_new;

        }
        return x;
    }
}

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
        res << (x[0] + x[1] - 6), (x[0] * x[0] - x[1]), (- x[0]), (- x[1]);
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


//    VectorXd L(VectorXd x, VectorXd lambda) {
//        VectorXd res(1);
//        res << f(x)[0] + lambda[1] * (x[0] * x[0] - x[1]) + lambda[0] * (x[0] + x[1] - 6) + lambda[2] * (6 - x[0]) +
//               lambda[3] * (6 - x[1]);
//        return res;
//    }
//
//    VectorXd nebula_x_L(VectorXd x, VectorXd lambda) {
//        VectorXd res(2);
//        res << 2 * (x[0] - 9.0 / 4.0) + lambda[0] + 2 * lambda[1] * x[0] - lambda[2], 2 * (x[1] - 2) + lambda[0] -
//                                                                                      lambda[1] - lambda[3];
//        return res;
//    }

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
            std::cout << "d:" << d << std::endl;
            std::cout << "alpha:" << alpha << std::endl;
            std::cout << "X:" << x << std::endl;
            printf("F:%f\n", f(x)[0]);

            x = x + alpha * d;
            lambda = lambda + alpha * (lambda_optimal - lambda);

        }
        return x;
    }
}

int main() {

# if 0
    test_wolfe_simplex_method();
    test_active_set_dual_method();
    test_primary_dual_interior_point_method();
#endif
# if 0
    Eigen::VectorXd x(2);
    x << 0, 0;
    Eigen::MatrixXd H(2, 2);
    H << 110, -6,
            -6, 18;
//    std::cout << Newton::Levenberg_Marquardt_solve(x, 4) << std::endl;
//    std::cout << Newton::BFGS_solve(x, H, 4) << std::endl;
    std::cout << Newton::DFP_solve(x, H.inverse(), 4) << std::endl;
#endif
# if 0
    Eigen::VectorXd x(2);
    x << 0, 0;
    Zigzag::line_search_solve(x, 10);
//    GoldenSearch::fibonacci_search_solve(0, 8);
#endif
    SQP::solve(10);
    return 0;

}


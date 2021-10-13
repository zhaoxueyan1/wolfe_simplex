
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

    return 0;
}


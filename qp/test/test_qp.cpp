//
// Created by zhaox on 2021/11/17.
//
#include <iostream>
#include "../wolfe_simplex/wolfe_simplex.h"
#include "../active_set_dual/QuadProg.h"
#include "../interior_point_qp/QuadraticProgramming.h"

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

# if 0
test_wolfe_simplex_method();
    test_active_set_dual_method();
    test_primary_dual_interior_point_method();
#endif
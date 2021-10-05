
#include <iostream>
#include "wolfe_simplex/wolfe_simplex.h"
#include "active_set_dual/QuadProg.h"

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
    SHG::Matdouble CE(2, 0, SHG::Vecdouble((std::vector<double>) {0,0}));
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
        printf("OK! The result is %f\n", f);
        for (int i = 0; i < x.size(); i++) {
            printf("X%d : %f ", i, x[i]);
        }
    } else {
        printf("Error %d\n", ok);
    }
}

int main() {
    test_wolfe_simplex_method();
    test_active_set_dual_method();
    test_primary_dual_interior_point_method();
    return 0;
}


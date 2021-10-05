
#include <iostream>
#include "wolfe_simplex.h"

int main() {
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
    return 0;
}


//
// Created by zhaox on 2021/12/1.
//

#include "simplex.h"

int main()
{
    int m, n;
    m = 14;//小等于约束的个数
    n = 6;//真实变量的个数,不包括松弛变量
    float st[10] = { 0 };
    for (int i = 0; i < 10; i++)
    {
        scanf("%f", &(st[i]));
    }
    for (int i = 0; i < 10; i++)
    {
        if (i % 2 == 0)
            //左区间
        {
            b_matrix[i] = -st[i];
            co_matrix[i][i / 2] = -1;
        }
        else
        {
            b_matrix[i] = st[i];
            co_matrix[i][i / 2] = 1;
        }
    }
    co_matrix[10][0] = 1, co_matrix[10][1] = -1, co_matrix[10][5] = 1, b_matrix[10] = 0;
    co_matrix[11][1] = 1, co_matrix[11][2] = -1, co_matrix[11][5] = 1, b_matrix[11] = 0;
    co_matrix[12][2] = 1, co_matrix[12][3] = -1, co_matrix[12][5] = 1, b_matrix[12] = 0;
    co_matrix[13][3] = 1, co_matrix[13][4] = -1, co_matrix[13][5] = 1, b_matrix[13] = 0;
    c_matrix[5] = -1;
    Simplex simplex(m, n);
    simplex.set_objective(c_matrix);
    simplex.set_co_matrix(co_matrix);
    simplex.set_bi_matrix(b_matrix);
    simplex.run();
    pair<vector<double>, double> rst = simplex.getans();
    printf("%0.2f\n", rst.second + 0.00005);
    for (int i = 0; i < rst.first.size()-1; i++)
    {
        printf("%0.2f", float(rst.first[i]+0.000005));
        if (i != (rst.first.size() - 2))
        {
            printf(" ");
        }
    }
    return 0;
}

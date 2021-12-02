//
// Created by zhaox on 2021/12/1.
//

#ifndef OPT_SIMPLEX_H
#define OPT_SIMPLEX_H


#include<iostream>
#include<algorithm>
#include<cstring>
#include<cstdio>
#include<cmath>
#include <vector>
using namespace std;
const double eps = 1e-20;
enum{ mxn = 50, mxm = 50 };//mxn是未知数的个数（包括松弛变量）mxm 约束个数
class Simplex
{
public:
    int n, m, t;
    double c[mxn] ;
    double a[mxm][mxn] ;
    int idx[mxn] , idy[mxn] ;
    int st[mxn] , top = 0;

    Simplex(int _m, int _n)
    {
        this->m = _m;
        this->n = _n;
        int i;
        for (i = 1; i <= n; i++)
        {
            idx[i] = i;//基变量
        }
        for (i = 1; i <= m; i++)
            //添加后的松弛变量
        {
            idy[i] = i + n;//非基变量
        }
    }
    void set_objective(double ci[mxn])
    {
        for (int i = 1; i <= n; i++)
        {
            a[0][i] = -ci[i - 1];//第0行作为目标函数
        }
    }
    void set_co_matrix(double co_matrix[mxm][mxn])
    {
        for (int i = 1; i <= m; i++)
        {
            for (int j = 1; j <= n; j++)
            {
                a[i][j] = co_matrix[i - 1][j - 1];
            }
        }
    }
    void set_bi_matrix(double b_matrix[mxm])
    {
        for (int i = 1; i <= m; i++)
        {
            a[i][0] = b_matrix[i-1];//第0列作为b矩阵
        }
        a[0][0] = 0;
    }
    int init_simplex()
    {
        while (1)
        {

            int i,x = 0, y = 0;
            for (i = 1; i <= m; i++)
            {
                if (a[i][0] < -eps && ((!x) || (rand() & 1)))
                    //如果有某个b[i] 即是约束小于0的
                {
                    x = i;
                }
            }
            if (!x)
                break;//没有小于0的
            for (i = 1; i <= n; i++)
            {
                if (a[x][i] < -eps && ((!y) || (rand() & 1)))
                    //从刚刚那一行bi小于0的那行,找到另外一个小于0的
                {
                    y = i;
                }
            }
            if (!y)
            {
                printf("Infeasible\n");
                return 0;
            }
            Pivot(x, y);//把第x行的第y列的元素作为主元 进行高斯消元
        }
        return 1;
    }
    void Pivot(int x, int y)
    {//用idy代换idx
        swap(idy[x], idx[y]);
        double tmp = a[x][y];
        a[x][y] = 1 / a[x][y];
        int i, j;
        top = 0;
        for (i = 0; i <= n; i++)
        {
            if (y != i)
                a[x][i] /= tmp;
        }

        for (i = 0; i <= n; i++)
        {
            if ((y != i) && fabs(a[x][i]) > eps)
            {
                st[++top] = i;
            }
        }
        for (i = 0; i <= m; i++)
        {
            if ((i == x) || (fabs(a[i][y]) < eps))
            {
                continue;
            }
            for (j = 1; j <= top; j++)
            {
                a[i][st[j]] -= a[x][st[j]] * a[i][y];
            }
            a[i][y] = -a[i][y] / tmp;
        }
        return;
    }
    int run()
    {
        int init=init_simplex();
        if (init == 0)
        {
            return init;//无解
        }
        int i, j;
        while (1){
            int x = 0, y = 0;
            double mn = 1e15;
            for (i = 1; i <= n; i++)
            {
                if (a[0][i] > eps)
                {
                    y = i;
                    break;
                }
            }
            if (!y)
                break;
            for (i = 1; i <= m; i++)
            {
                if (a[i][y] > eps && (a[i][0] / a[i][y] < mn))
                {
                    mn = a[i][0] / a[i][y];
                    x = i;
                }
            }
            if (!x)
            {
                printf("Unbounded\n");
                return -1;//无界
            }
            Pivot(x, y);
        }
        return 1;//有解
    }
    pair<vector<double>, double> getans()
    {
        vector<double> x;
        double z;
        int i;
        z = a[0][0];
        for (i = 1; i <= n; i++)
        {
            a[0][i] = 0;
        }
        for (i = 1; i <= m; i++)
        {
            if (idy[i] <= n)a[0][idy[i]] = a[i][0];
        }
        for (i = 1; i <= n; i++)
        {
            x.push_back(a[0][i]);
        }
        return pair< vector<double>, double>(x, z);
    }
};
double b_matrix[mxm] = { 0 };
double co_matrix[mxm][mxn] = { 0 };
double c_matrix[mxn] = { 0 };



#endif //OPT_SIMPLEX_H

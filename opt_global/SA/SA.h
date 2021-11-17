//
// Created by zhaox on 2021/11/17.
//

#ifndef OPT_SA_H
#define OPT_SA_H
#include <Eigen/Dense>
#include <cmath>

namespace SA {
    using namespace Eigen;
    double e=1e-5,at=0.98,T=1.0;
    double PI = 3.1415926535;
    double f(double x) {
        double res = 0;
        res = -exp(-2* log(2)*((x-0.008)/0.854)*pow(sin(5*PI*(pow(x,0.75)-0.05)),6));
        return res;
    }
    //e表示终止温度  at表示温度的变化率  T是初始温度
//    int L = 200000; //L为最大迭代次数
    void Simulated_Annealing(int itr) {
        int L = itr;
        //模拟退火
        double x = 0;
        double pre_f = f(0);
        while (L--) {  //最多迭代L次
            double new_x = x + (1.*(rand()%INT16_MAX)-INT16_MAX/2) / (INT16_MAX);
            if(x<0 | x > 1){
                continue;
            }
            double df = f(new_x) - pre_f;  //计算代价
            double sj = rand() % 10000;     //sj为0~1之间的随机数
            printf("%f\t%f\t%f\n",x,pre_f,exp(-df / T));
            sj /= 10000;
            if (df < 0) {  //如果结果更优，直接接受
                x = new_x;
                pre_f += df;
            } else if (exp(-df / T) > sj) {  //否则以一定概率接受
                x = new_x;
                pre_f += df;
            }
            T *= at;  //降温
            if (T < e) break;  //T降到终止温度后退出
        }
    }
}


#endif //OPT_SA_H

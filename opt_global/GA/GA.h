//
// Created by zhaox on 2021/11/17.
//

#ifndef OPT_GA_H
#define OPT_GA_H

//功能：求一个多元函数的最大值（这里是求二元函数的最大值：f(x1,x2) = 21.5+x1*sin(4pi*x1)+x2*sin(20pi*x2)）
#pragma once//保证文件只被编译一次

#include<random>
#include<vector>
#include<iostream>
#include<cmath>
#include<ctime>
#include <cstdlib>
#include <bitset>
#include<iomanip>

using namespace std;

const double PI = 3.141592653589793;//定义一个不可改变的常量值PI
const int Po_Size = 300;//种群规模
const int Ev_Algebra = 5000;//进化代数
const double Ov_Probability = 0.850; //交叉概率,交叉概率用于判断两两个体是否需要交叉
const double Va_Probability = 0.50;//变异概率,变异概率用于判断任一个体是否需要变异
const int De_Variable = 2;//函数自变量的个数,如果要进行多元函数的最值求解，直接修改自变量数个De_Variable即可
const int length1 = 24;//精确到6位小数，用24位二进制编码
const int length2 = 23;//精确到6位小数，用23位二进制编码

class X_Range //自变量取值范围类，适用于多个变量
{
private:
    double Upper;//变量的上界取值
    double Lower;//变量的下界取值
public:
    X_Range(double m_Upper, double m_Lower);//构造函数
    double GetUpper() const;//获取变量上限
    double GetLower() const;//获取变量下限
};

class Individual //定义个体类
{
private:
    double Variable[De_Variable];//存放变量x1,x2,x3........
    double Fitness;//适应值
    double ReFitness;//适应值概率
    double SumFitness;//累加概率，为轮盘转做准备
public:
    Individual() {}//默认构造函数
    Individual(double *m_Variable);//构造函数
    double *GetVariable();//获取变量值
    void ChaFitness(const double m_fitness);//修改适应值
    void ChaReFitness(const double m_ReFitness);//修改适应值概率
    void ChaSumFitness(const double m_SumFitness);//修改累加概率
    double GetFitness() const;//获取适应值
    double GetReFitness() const;//获取适应值概率
    double GetSumFitness() const;//获取累加概率
};

void Initialize();//随机初始化种群，得到第一代个体
void CaculaFitness();//计算个体的适应值
void CaculaReFitness();//计算个体的适应值概率
void CalculaSumFitness();//计算累加个体概率
void seclect();//种群选择
double Scand();//随机产生0到49的随机整数
void crossing();//杂交
void variating();//变异
void genetic_algorithm();//遗传算法

#endif //OPT_GA_H

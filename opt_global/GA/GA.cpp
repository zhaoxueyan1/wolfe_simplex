//
// Created by zhaox on 2021/11/17.
//

#include "GA.h"

//自变量取值范围向量和种群向量定义
const X_Range Range[De_Variable] = {X_Range(-6.0, 6.0), X_Range(-6.0, 6.0)};//自变量（或者基因）x1,x2的取值范围
vector<Individual> nowpopulation;//P(t)种群
vector<Individual> midpopulation;//中间种群，存放轮盘选择后的优秀个体
vector<Individual> nextpopulation;//P(t+1)种群
//X_Range类实现
X_Range::X_Range(double m_Lower, double m_Upper) : Lower(m_Lower), Upper(m_Upper) {}//X_Range类构造函数实现
double X_Range::GetUpper() const//获取变量上限
{
    return Upper;
}

double X_Range::GetLower() const//获取变量下限
{
    return Lower;
}

//Individual类实现
Individual::Individual(double *m_Variable)//构造函数
{
    for (int i = 0; i < De_Variable; i++)//用for循环自变量逐个赋值
    {
        if (m_Variable[i] >= Range[i].GetLower() && m_Variable[i] <= Range[i].GetUpper())//这里要进行自变量取值范围判断
        {
            Variable[i] = m_Variable[i];//自变量赋值
        } else//不满足要求则发出出错警告并返回
        {
            cerr << "自变量取值不满足要求" << endl;
            exit(1);//停止程序，我会以随机函数的方式生成自变量的值(基因值)，这里说明基因值不在规定范围内
        }
    }
    //初始化时默认适应值等值为0
    this->Fitness = 0;
    this->ReFitness = 0;
    this->SumFitness = 0;
}

double *Individual::GetVariable()//获取基因值
{
    return Variable;
}

double Individual::GetFitness() const//获取适应值
{
    return Fitness;
}

double Individual::GetReFitness() const //获取适应值概率
{
    return ReFitness;
}

double Individual::GetSumFitness() const//获取累加概率
{
    return SumFitness;
}

void Individual::ChaFitness(const double m_fitness)//修改适应值
{
    this->Fitness = m_fitness;
}

void Individual::ChaReFitness(const double m_ReFitness)//修改适应值概率
{
    this->ReFitness = m_ReFitness;
}

void Individual::ChaSumFitness(const double m_SumFitness)//修改累加概率
{
    this->SumFitness = m_SumFitness;
}

//遗传算法的准备工作
void Initialize()//随机初始化种群，得到第一代种群
{
//产生指定范围的随机变量（基因）
    double X[Po_Size][De_Variable];//为了使程序可以满足多元函数最值的计算，用矩阵保存产生的随机数变量值
    for (int j = 0; j < De_Variable; j++) {
        default_random_engine e(time(0));//引擎，生成随机序列
        uniform_real_distribution<double> u(Range[j].GetLower(), Range[j].GetUpper());//分布
        for (int i = 0; i < Po_Size; i++)//先按列存储随机数
        {
            X[i][j] = u(e);//循环结束时，所有随机值就保存在X矩阵中
        }
    }
//生成对象（染色体）并加入到初始种群中
    for (int i = 0; i < Po_Size; i++) {
        double variable[De_Variable];
        for (int j = 0; j < De_Variable; j++) {
            variable[j] = X[i][j];//按行保存
        }
        Individual Indivi(variable);//生成一个对象（染色体）
        nowpopulation.push_back(Indivi);//加入到种群population中
    }
}

void CaculaFitness()//计算个体的适应值
{
    //f(x1,x2) = 21.5+x1*sin(4pi*x1)+x2*sin(20pi*x2)）为适应度计算函数
    double fitness = 0;//临时适应值
    double x[De_Variable];//临时存储自变量（基因）
    for (int i = 0; i < Po_Size; i++) {
        for (int j = 0; j < De_Variable; j++)
            x[j] = nowpopulation.at(i).GetVariable()[j];//这样更直观
        fitness = (2186.0 - pow(x[0] * x[0] + x[1] - 11, 2) \
 - pow(x[0] + x[1] * x[1] - 7, 2)) / 2186.0;//适应度计算
        nowpopulation.at(i).ChaFitness(fitness);//修改当前染色体的适应值
    }
}

void CaculaReFitness()//计算适应值概率
{
    double sum = 0;//适应值累加器
    double temp = 0;
    for (int i = 0; i < Po_Size; i++)//计算出适应值之和
    {
        sum += nowpopulation.at(i).GetFitness();
    }
    for (int j = 0; j < Po_Size; j++) {
        temp = nowpopulation.at(j).GetFitness() / sum;//计算概率
        nowpopulation.at(j).ChaReFitness(temp);//修改个体的适应度概率
    }
}

void CalculaSumFitness()//计算累加个体概率
{
    double summation = 0;//累加器
    for (int k = 0; k < Po_Size; k++) {
        summation += nowpopulation.at(k).GetReFitness();
        nowpopulation.at(k).ChaSumFitness(summation);//当前累加结果赋值
    }
}

void seclect() //种群选择
{
    //随机生生成0到1的小数
    double array[Po_Size];//随机数保存变量
    default_random_engine e(time(0));//引擎，生成随机序列
    uniform_real_distribution<double> u(0.0, 1.0); //分布
    for (int i = 0; i < Po_Size; i++)
        array[i] = u(e);
    //轮盘进行选择
    for (int j = 0; j < Po_Size; j++) {
        for (int i = 1; i < Po_Size; i++) {
            if (array[j] < nowpopulation[i - 1].GetSumFitness()) {
                midpopulation.push_back(nowpopulation.at(i - 1));//加入到中间种群
            }
            if (array[j] >= nowpopulation.at(i - 1).GetSumFitness() &&
                array[j] <= nowpopulation.at(i).GetSumFitness()) {
                midpopulation.push_back(nowpopulation.at(i));//加入到中间种群
            }
        }
    }
    nowpopulation.clear();//清空nowpopulation
}

double Scand() //随机产生0到1的小数
{
    int N = rand() % 999;
    return double(N) / 1000.0;;//随机产生0到1的小数
}

void crossing()//杂交
{
    int num = 0;//记录次数
    double corss = 0.0;//保存随机产生的概率值
    srand((unsigned) time(NULL));//根据系统时间设置随机数种子,设置一次随机种子就行
    double array1[De_Variable], array2[De_Variable];//临时存储父亲和母亲的变量值
    while (num < Po_Size - 1)//个体1与个体2杂交，个体3与个体4杂交......个体i和个体i+1杂交
    {
        //判断双亲是否需要杂交，随机生成一个0到1的小数，如果这个数大于杂交概率，则放弃杂交，直接遗传给下一代，否则，对父母体进行杂交
        corss = Scand();
        if (corss <= Ov_Probability)//如果corss小于等于杂交概率Ov_Probability就进行单点杂交
        {
            //首先寻找对应下标的个体并且保存
            for (int i = 0; i < De_Variable; i++) {
                array1[i] = midpopulation.at(num).GetVariable()[i];//父亲的自变量
                array2[i] = midpopulation.at(num + 1).GetVariable()[i];//母亲自变量
            }
            int localx1, localx2;//记录基因交叉点的位置
            int corssx1[length1], corssx2[length2];//作为交换基因的数组
            double newx1[2], newx2[2];//分别用来保存基因交换后所对应自变量值
            bool p1 = true, p2 = true;
            //然后对双亲变量进行编码并且进行单点杂交
            for (int j = 0; j < De_Variable; j++)//array1的x1编码之后和array2的x1编码后进行单点杂交，以此类推
            {
                if (j == 0)//x1进行编码并且杂交
                {
                    bitset <length1> array1b1((array1[j] + 3.0) * pow(10, 6));//加上3.0形成一个unsigaed数之后在进行母体1的x1编码
                    bitset <length1> array2b1((array2[j] + 3.0) * pow(10, 6));//加上3.0形成一个unsigaed数之后在进行母体2的x1编码
                    //现在随机生成0到length1-1的数，确定交叉点的位置
                    localx1 = rand() % length1;
                    //现在进行单点交叉，交换双亲localx1后面的基因
                    for (int i = 0; i < localx1; i++)
                        corssx1[i] = array1b1[i];
                    for (int k = 0; k < localx1; k++)
                        array1b1[k] = array2b1[k];
                    for (int s = 0; s < localx1; s++)
                        array2b1[s] = corssx1[s];
                    //新值保存在newx1数组中，x1基因完成单点杂交操作
                    newx1[0] = double(array1b1.to_ullong()) / pow(10, 6) - 3.0;
                    newx2[0] = double(array2b1.to_ullong()) / pow(10, 6) - 3.0;
                    //对新产生的值进行判断，判断是否超出范围，如果超出范围则不杂交
                    if (newx1[0] < Range[0].GetLower() || newx1[0] > Range[0].GetUpper() ||
                        newx2[0] < Range[0].GetLower() || newx2[0] > Range[0].GetUpper()) {
                        p1 = false;
                        break;
                    }
                }
                if (j == 1)//x2进行编码并且杂交
                {
                    bitset <length2> array1b2((array1[j]) * pow(10, 6));//母体1的x2编码
                    bitset <length2> array2b2((array2[j]) * pow(10, 6));//母体2的x2编码
                    //现在随机生成0到length2-1的数，确定交叉点的位置
                    localx2 = rand() % length2;
                    //现在进行单点交叉，交换双亲localx2后面的基因
                    for (int i = 0; i < localx2; i++)
                        corssx2[i] = array1b2[i];
                    for (int k = 0; k < localx2; k++)
                        array1b2[k] = array2b2[k];
                    for (int s = 0; s < localx2; s++)
                        array2b2[s] = corssx2[s];
                    //新值保存在newx2数组中，x2基因完成单点杂交操作
                    newx1[1] = double(array1b2.to_ullong()) / pow(10, 6);
                    newx2[1] = double(array2b2.to_ullong()) / pow(10, 6);
                    //对新产生的值进行判断，判断是否超出范围，如果超出范围则不杂交
                    if (newx1[1] < Range[1].GetLower() || newx1[1] > Range[1].GetUpper() ||
                        newx2[1] < Range[1].GetLower() || newx2[1] > Range[1].GetUpper()) {
                        p2 = false;
                        break;
                    }
                }
            }
            if (p1 == true && p2 == true) {
                Individual newchiled1(newx1);
                Individual newchiled2(newx2);
                nextpopulation.push_back(newchiled1);
                nextpopulation.push_back(newchiled2);
            } else//将原来的个体遗传给下一代
            {
                nextpopulation.push_back(midpopulation.at(num));
                nextpopulation.push_back(midpopulation.at(num + 1));
            }
        } else//否则直接遗传给下一代nextpopulation
        {
            nextpopulation.push_back(midpopulation.at(num));//生成一个新的个体并且加入到nextpopulation中
            nextpopulation.push_back(midpopulation.at(num + 1));
        }
        num += 2;
    }
    midpopulation.clear();//清空midpopulation
}

void variating()//变异
{
    int num = 0;
    while (num < Po_Size) {
        double variation = Scand();//随机产生一个0到1的小数，用于判断是否进行变异
        if (variation <= Va_Probability)//如果variation小于变异系数，则需要进行变异
        {
            double x[2];
            bool p = true;
            int x1local, x2local;
            x[0] = nextpopulation.at(num).GetVariable()[0];
            x[1] = nextpopulation.at(num).GetVariable()[1];
            bitset <length1> array1((x[0] + 3.0) * pow(10, 6));//x1编码
            bitset <length2> array2(x[1] * pow(10, 6));//x2编码
            x1local = rand() % length1;//array1该位取反
            x2local = rand() % length2;//array2该位取反
            array1.flip(x1local);//改变array1 x1local位的状态
            array2.flip(x2local);//改变array2 x2local位的状态
            x[0] = double(array1.to_ullong()) / pow(10, 6) - 3.0;
            x[1] = double(array2.to_ullong()) / pow(10, 6);
            //判断是否符合条件
            if (x[0] < Range[0].GetLower() || x[0] > Range[0].GetUpper() || x[1] < Range[1].GetLower() ||
                x[1] > Range[1].GetUpper())
                p = false;
            if (!p)
                nowpopulation.push_back(nextpopulation.at(num));
            if (p) {
                Individual newchiled(x);
                nowpopulation.push_back(newchiled);
            }
        } else
            nowpopulation.push_back(nextpopulation.at(num));
        num++;
    }
    nextpopulation.clear();//清空nextpopulation
}

void genetic_algorithm() {
    Initialize();//初始化种群,随机生成第一代个体
    //进化500代
    for (int i = 0; i < Ev_Algebra; i++) {
        CaculaFitness();//适应度计算
        CaculaReFitness();//适应度概率计算
        CalculaSumFitness();//计算累加个体概率
        seclect();//选择
        crossing();//杂交
        variating();//变异
    }
    CaculaFitness();//适应度计算
    double maxfitness = nowpopulation.at(0).GetFitness();
    int maxid = 0;
    int k;
    for (k = 0; k < Po_Size; k++) {
        if (maxfitness < nowpopulation.at(k).GetFitness()) {
            maxfitness = nowpopulation.at(k).GetFitness();
            maxid = k;
        }
    }
    //进化500代之后输出
    cout << "x1" << setw(10) << "x2" << setw(15) << "Fitness" << endl;
    for (int j = 0; j < Po_Size; j++)
        cout << nowpopulation.at(j).GetVariable()[0] << setw(10) << nowpopulation.at(j).GetVariable()[1] << setw(10)
             << nowpopulation.at(j).GetFitness() << endl;
    cout << "x1=" << nowpopulation.at(maxid).GetVariable()[0] << " ，" << "x2="
         << nowpopulation.at(maxid).GetVariable()[1] << "时取得最大值：" << maxfitness << endl;
}
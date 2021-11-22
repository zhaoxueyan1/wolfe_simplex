//
// Created by zhaox on 2021/11/22.
//

#include "bb.h"
#include<iostream>
#include <queue>
#include <cmath>
using namespace std;

struct node{
    int depth;
    double now_val;
    double pre_y;
    vector<int> u;

};

double next_y(double now_y,double now_u){
    return 0.9*now_y+0.1*now_u;
}
int du[]={0,1};

double BB(){
    queue<node>Q;
    node s;
    s.depth = 0;
    s.now_val = 0;
    s.pre_y = 0.5;
    int max_depth = 3;
    double MIN_VAL = 1e12;
    vector<int> res;
    Q.push(s);
    while(!Q.empty()){
        node t = Q.front();Q.pop();
        if (t.depth == max_depth){
            if(t.now_val <= MIN_VAL){
                MIN_VAL = t.now_val;
                res = t.u;
            }
            continue;
        }

        for (int i = 0; i < 2; ++i) {
            int u = du[i];
            vector<int> now_u = t.u;
            now_u.push_back(u);
            int now_depth = t.depth +1;
            double now_y = next_y(t.pre_y,u);
            double now_val = t.now_val + pow(now_y -1,2.0) + 0.01*u;
            double all_val = now_val + (max_depth -now_depth)*(1 + 0.01*1);
            if(all_val <= MIN_VAL){
                MIN_VAL = all_val;
                printf("%d\t%d\t%f\t\n",now_depth,u,pow(now_y -1,2.0) + 0.01*u);
                Q.push({now_depth,now_val,now_y,now_u});
            }
        }
    }
    for (auto i:res) {
        printf("%d ",i);
    }
    putchar(10);
    return MIN_VAL;
}
double dp(){


}
int main(){
    printf("Min val %f\n",BB());

    return 0;
}
//
// Created by zhaox on 2021/10/5.
//

#ifndef OPT_QUADRATICPROGRAMMING_H
#define OPT_QUADRATICPROGRAMMING_H

#pragma once
#include <Eigen/Dense>
#include <vector>
//#include "mathOperator.h"

class QuadraticProgramming
{
public:
    enum Algorithm
    {
        PrimalDualInteriorPoint
    };
    enum ExitFlag
    {
        SolutionFound,
        InfeasibleProblem,
        ExceedMaxIternation,
        WrongDimension
    };

    QuadraticProgramming(const Algorithm& algo = (Algorithm::PrimalDualInteriorPoint));
    QuadraticProgramming(const Eigen::MatrixXd& P, const Eigen::VectorXd& q, const Eigen::MatrixXd& G = (Eigen::MatrixXd()), const Eigen::VectorXd& h = (Eigen::VectorXd()),
                         const Eigen::MatrixXd& A = (Eigen::MatrixXd()), const Eigen::VectorXd& b = (Eigen::VectorXd()), const Algorithm& algo = (Algorithm::PrimalDualInteriorPoint));

    Eigen::VectorXd solve(const Eigen::VectorXd& x0);
    Eigen::VectorXd PrimalDualInteriorPointMethod(const Eigen::VectorXd& x0);
    Eigen::VectorXd projectToFeasibleSpace(const Eigen::VectorXd& x0);
    Algorithm _algorithm;
    ExitFlag _exit;

    // objective function: 0.5*x'*P*x + q'*x
    Eigen::MatrixXd _P;
    Eigen::VectorXd _q;
    // inequality constraint: G*x <= h
    Eigen::MatrixXd _G;
    Eigen::VectorXd _h;
    // equality constraint: A*x = b ?
    Eigen::MatrixXd _A;
    Eigen::VectorXd _b;

    Eigen::VectorXd _xf;		// primal variable
    Eigen::VectorXd _lambda;	// dual variable for ineq
    Eigen::VectorXd _nu;		// dual variable for eq

    unsigned int _Iter;
    unsigned int _maxIteration;
    int _iterProj;
    int _maxIterationProj;
    double _tolFea;
    double _eps;
    bool _Display;
    // line search parameters
    double _alpha;
    double _beta;
    // primal dual iteration parameter
    double _mu;
};
#endif //OPT_QUADRATICPROGRAMMING_H

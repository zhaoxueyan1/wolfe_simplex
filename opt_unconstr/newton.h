//
// Created by zhaox on 2021/10/10.
//

#ifndef OPT_NEWTOWN_H
#define OPT_NEWTOWN_H

#include <Eigen/Dense>

/*
 * f(x) = 1/2x^T *H *x + p^T * x
 * \nebula f(x) = H x + p
 * \Hession f = H
 * */
class Newton {
public:
    Newton(const Eigen::MatrixXd &H, const Eigen::MatrixXd &p)
            : _H(H), _p(p) {}

    Eigen::VectorXd NewTonSolve() {
        Eigen::VectorXd res;
        for (int i = 0; i < _p.count(); ++i) {
            res << 0;
        }
        int itr = _itr;
        while (--itr) {
            res = res - _H.inverse()* (_H* res + _p);
        }
        return res;
    }
    Eigen::VectorXd LevenbergSolve() {
        Eigen::VectorXd res;
        for (int i = 0; i < _p.count(); ++i) {
            res << 0;
        }
        int itr = _itr;
        while (--itr) {
            res = res - _H.inverse()* (_H* res + _p);
        }
        return res;
    }
    Eigen::VectorXd BroydenSolve() {
        Eigen::VectorXd res;
        for (int i = 0; i < _p.count(); ++i) {
            res << 0;
        }
        int itr = _itr;
        while (--itr) {
            res = res - _H.inverse()* (_H* res + _p);
        }
        return res;
    }
    Eigen::VectorXd DavidonSolve() {
        Eigen::VectorXd res;
        for (int i = 0; i < _p.count(); ++i) {
            res << 0;
        }
        int itr = _itr;
        while (--itr) {
            res = res - _H.inverse()* (_H* res + _p);
        }
        return res;
    }
private:
    Eigen::MatrixXd _H;
    Eigen::VectorXd _p;
    static const int _itr = 5;

};


#endif //OPT_NEWTOWN_H

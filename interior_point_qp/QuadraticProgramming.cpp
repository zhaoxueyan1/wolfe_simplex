#include <iostream>
#include "QuadraticProgramming.h"

using namespace std;

QuadraticProgramming::QuadraticProgramming(const QuadraticProgramming::Algorithm& algo) : QuadraticProgramming(Eigen::MatrixXd(), Eigen::VectorXd(), Eigen::MatrixXd(), Eigen::VectorXd(), Eigen::MatrixXd(), Eigen::VectorXd(), algo) {}

QuadraticProgramming::QuadraticProgramming(const Eigen::MatrixXd& P, const Eigen::VectorXd& q, const Eigen::MatrixXd& G, const Eigen::VectorXd& h, const Eigen::MatrixXd& A, const Eigen::VectorXd& b, const Algorithm& algo) :
        _P(P), _q(q), _G(G), _h(h), _A(A), _b(b), _Display(false)
{
    if (algo == QuadraticProgramming::Algorithm::PrimalDualInteriorPoint)
    {
        _algorithm = QuadraticProgramming::Algorithm::PrimalDualInteriorPoint;
        _maxIteration = 300;
        _maxIterationProj = 200;
        _tolFea = 5e-5;
        _eps = 5e-5;
        _mu = 10; // >=10
        _beta = 0.5; // 0.3~0.8
        _alpha = 0.01; // 0.01~0.1
        _Display = true;
    }
}

Eigen::VectorXd QuadraticProgramming::solve(const Eigen::VectorXd& x0)
{
    if (_algorithm == QuadraticProgramming::PrimalDualInteriorPoint)
    {
        return PrimalDualInteriorPointMethod(x0);
    }
    return x0;
}

Eigen::VectorXd QuadraticProgramming::PrimalDualInteriorPointMethod(const Eigen::VectorXd& x0)
{

    if (_G.rows() != _h.size() || _A.rows() != _b.size())
    {
        _exit = ExitFlag::WrongDimension;
        return x0;
    }
    if (_G.rows() > 0)
    {
        if (x0.size() != _G.cols())
        {
            _exit = ExitFlag::WrongDimension;
            return x0;
        }
    }
    if (_A.rows() > 0)
    {
        if (x0.size() != _A.cols())
        {
            _exit = ExitFlag::WrongDimension;
            return x0;
        }
    }
    Eigen::VectorXd residual(x0.size() + _A.rows());
    Eigen::VectorXd rdual(x0.size());
    Eigen::VectorXd rpri(_A.rows());
    Eigen::VectorXd rcent(_G.rows());
    Eigen::VectorXd f_i(_G.rows());
    Eigen::VectorXd dy(x0.size() + _A.rows());
    Eigen::VectorXd dlambda(_G.rows());
    Eigen::MatrixXd Mtx_pd(x0.size() + _A.rows(), x0.size() + _A.rows());
    Eigen::VectorXd x_cur;
    Eigen::VectorXd lambda_cur;
    Eigen::VectorXd nu_cur;
    Eigen::VectorXd tmp;
    Mtx_pd.setZero();
    if (_A.rows() > 0)
    {
        Mtx_pd.block(x0.size(), 0, _A.rows(), x0.size()) = _A;
        Mtx_pd.block(0, x0.size(), x0.size(), _A.rows()) = _A.transpose();
    }
    bool isFeasible;
    double tinv;
    double t;
    double s;
    double etha;
    double r_next;
    double r_cur;

    // initialize
    if ((_G.rows() > 0 && (_G*x0 - _h).maxCoeff() > - _eps) || (_A.rows() > 0 && (_A*x0 - _b).squaredNorm() > 10.0*_tolFea))
    {
        _xf = projectToFeasibleSpace(x0);
        //cout << _xf.transpose() << endl;
        if (_iterProj >= _maxIterationProj)
        {
            _exit = ExitFlag::InfeasibleProblem;

            //cout << "G=" << endl << _G << endl;
            //cout << "h=" << endl << _h << endl;
            //cout << "A=" << endl << _A.transpose() << endl;
            //cout << "b=" << endl << _b << endl;
            return _xf;
        }
        if ((_A.rows() > 0) && (_A*_xf - _b).norm() > 10.0*_tolFea)
        {
            _exit = ExitFlag::InfeasibleProblem;
            return _xf;
        }
    }
    else
        _xf = x0;
    _lambda = Eigen::VectorXd::Ones(_G.rows());
    _nu = Eigen::VectorXd::Zero(_A.rows());
    rdual.setOnes();
    rpri.setOnes();
    if (_G.rows() > 0)
        etha = 1.0;
    else
        etha = 0.0;
    _Iter = 0;

    vector<double> f_i_vec(_G.rows());
    while (_Iter < _maxIteration && !(rdual.squaredNorm() < _tolFea && rpri.squaredNorm() < _tolFea && etha < _eps))
    {
        _Iter++;
        if (_G.rows() > 0)
        {
            f_i = _G*_xf - _h;
            // for debugging
            for (int i = 0; i < _G.rows(); i++)
                f_i_vec[i] = f_i(i);
            etha = -f_i.transpose()*_lambda;
            t = _mu*_G.rows() / etha;
            tinv = 1.0 / t;
        }
        // find line search direction
        rdual = _P*_xf + _q;
        if (_A.rows() > 0)
        {
            rdual += _A.transpose()*_nu;
            rpri = _A*_xf - _b;
        }
        residual.segment(0, x0.size()) = rdual;
        if (_G.rows() > 0)
            rdual += _G.transpose()*_lambda;
        for (int i = 0; i < _G.rows(); i++)
        {
            residual.segment(0, x0.size()) -= tinv*_G.row(i).transpose() / f_i(i);
            //printf("%f\t%f\n", residual(0), f_i(i));
        }

        residual.segment(x0.size(), _A.rows()) = rpri;
        residual *= -1.0;
        Mtx_pd.block(0, 0, x0.size(), x0.size()) = _P;
        for (int i = 0; i < _G.rows(); i++)
        {
            Mtx_pd.block(0, 0, x0.size(), x0.size()) -= _lambda(i) / f_i(i) * _G.row(i).transpose()*_G.row(i);
        }
        //cout << "Matrix: " << endl << Mtx_pd << endl;
        dy = Mtx_pd.ldlt().solve(residual);
        //cout << (residual - Mtx_pd*dy).transpose() << endl;
        //cout << "dy: " << dy.transpose() << endl;
        if (_G.rows() > 0)
            dlambda = _G*dy.segment(0, x0.size());
        for (int i = 0; i < _G.rows(); i++)
        {
            rcent(i) = -_lambda(i)*f_i(i) - tinv;
            dlambda(i) *= -_lambda(i);
            dlambda(i) += rcent(i);
            dlambda(i) /= f_i(i);
        }

        // line search
        s = 1.0;
        //cout << "l: " << _lambda.transpose() << endl;
        //cout << "dl: " << dlambda.transpose() << endl;
        for (int i = 0; i < _G.rows(); i++)
        {
            if (dlambda(i) < 0)
                s = min(s, -_lambda(i) / dlambda(i));
        }
        s = 0.99*s;

        r_cur = sqrt(rdual.squaredNorm() + rpri.squaredNorm() + rcent.squaredNorm());
        x_cur = _xf;
        lambda_cur = _lambda;
        nu_cur = _nu;

        // reduce s to stay within feasible region
        int test = 0;
        if (_G.rows() > 0)
        {
            while (test < 300)
            {
                test++;
                isFeasible = false;
                _xf = x_cur + s*dy.segment(0, x0.size());
                tmp = _G*_xf - _h;
                //cout << tmp.transpose() << endl;
                if (tmp.maxCoeff() < 0.0)
                    isFeasible = true;
                if (isFeasible)
                    break;
                s *= _beta;
            }
            if (test == 300)
            {
                _exit = ExitFlag::ExceedMaxIternation;
                return _xf;
            }
        }
        int test2 = 0;
        while (1)
        {
            test2++;
            _xf = x_cur + s*dy.segment(0, x0.size());
            _lambda = lambda_cur + s*dlambda;
            _nu = nu_cur + s*dy.segment(x0.size(), _A.rows());
            rdual = _P*_xf + _q;
            if (_A.rows() > 0)
            {
                rdual += _A.transpose()*_nu;
                rpri = _A*_xf - _b;
            }
            if (_G.rows() > 0)
                rdual += _G.transpose()*_lambda;
            for (int i = 0; i < _G.rows(); i++)
                rcent(i) = -_lambda(i)*f_i(i) - tinv;
            r_next = sqrt(rdual.squaredNorm() + rpri.squaredNorm() + rcent.squaredNorm());
            if (r_next <= (1 - _alpha*s)*r_cur)
                break;
            s *= _beta;
            if (test2 == 300)
            {
                _exit = ExitFlag::ExceedMaxIternation;
                return _xf;
            }
        }
        //if (_Display && _Iter % 10 == 0)
        printf("iter = %u\n",_Iter);
        std::cout<<_xf<<endl;
//        printf("X = %.3f, iter = %u\n", 0.5*_xf.transpose()*_P*_xf + _xf.transpose()*_q, _Iter);
    }
    if (_Iter >= _maxIteration)
        _exit = ExitFlag::ExceedMaxIternation;
    else
        _exit = ExitFlag::SolutionFound;
    return _xf;
}

Eigen::VectorXd QuadraticProgramming::projectToFeasibleSpace(const Eigen::VectorXd & x0)
{
    Eigen::MatrixXd Aug(_A.rows() + _G.rows(), x0.size() + _G.rows());
    Eigen::VectorXd vAug(_A.rows() + _G.rows());
    Aug.setZero();
    if (_A.rows() > 0)
    {
        Aug.block(0, 0, _A.rows(), x0.size()) = _A;
        vAug.segment(0, _A.rows()) = _b;
    }
    if (_G.rows() > 0)
    {
        Aug.block(_A.rows(), 0, _G.rows(), x0.size()) = _G;
        vAug.segment(_A.rows(), _G.rows()) = _h;
    }

    for (int i = 0; i < _G.rows(); i++)
        Aug(_A.rows() + i, x0.size() + i) = 1.0;

    Eigen::VectorXd xNs = Aug.transpose()*(Aug*Aug.transpose()).ldlt().solve(vAug);
    //cout << (Aug*xNs - vAug).transpose() << endl;
    Eigen::MatrixXd N = Eigen::MatrixXd::Identity(x0.size() + _G.rows(), x0.size() + _G.rows()) - Aug.transpose()*(Aug*Aug.transpose()).inverse()*Aug;
    xNs += N.block(0, 0, xNs.size(), x0.size())*x0;
    _iterProj = 0;
    //cout << (Aug*xNs - vAug).transpose() << endl;
    //cout << "checkMatrix: " << endl;
    //cout << Aug*Aug.transpose() << endl;
    //cout << "checkMatrix: " << endl;
    //cout << (Aug*Aug.transpose()).inverse() << endl;
    if (_G.rows() > 0)
    {
        bool isFeasible = false;
        double step_size = 0.1;
        Eigen::VectorXd updateDir(xNs.size());

        while (!isFeasible && _iterProj < _maxIterationProj)
        {
            _iterProj++;
            updateDir.setZero();
            for (int i = 0; i < _G.rows(); i++)
            {
                if (xNs(x0.size() + i) < _eps)
                    updateDir(x0.size() + i) = -step_size*min(xNs(x0.size() + i), 0.0) + step_size;
            }
            xNs += N*updateDir;
            if (xNs.tail(_G.rows()).minCoeff() > _eps)
                isFeasible = true;
        }
    }
    // for debugging
    vector<double> checkproj(_G.rows());
    Eigen::VectorXd check = (_G*xNs.head(x0.size()) - _h);
    for (int i = 0; i < _G.rows(); i++)
        checkproj[i] = check(i);
    vector<double> checkproj2(_A.rows());

    Eigen::VectorXd check2 = Eigen::VectorXd();
    if (_A.rows() > 0)
    {
        check2 = (_A*xNs.head(x0.size()) - _b);
        if (check2.norm() > _eps)
            printf("singular problem\n");
    }

    for (int i = 0; i < _A.rows(); i++)
        checkproj2[i] = check2(i);

    vector<double> s(_G.rows());
    for (unsigned int i = 0; i < s.size(); i++)
        s[i] = xNs(x0.size() + i);
    //cout << "checkproj: " << (_G*xNs.head(x0.size()) - _h).transpose() << endl;

    return xNs.head(x0.size());
}
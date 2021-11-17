#pragma once

#define EPS_DLSI			1.0e-2

#include <cmath>
#include <vector>
#include <iostream>

#define lambdamax_DLSI		0.2
//#define std::min(a, b)  (((a) < (b)) ? (a) : (b))


static Eigen::MatrixXd pinv(const Eigen::MatrixXd& _matrix)
{
    Eigen::MatrixXd singInvMatrix(_matrix.cols(), _matrix.rows());
    singInvMatrix.setZero();
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(_matrix, Eigen::ComputeFullU | Eigen::ComputeFullV);
    for (int i = 0; i < std::min(_matrix.cols(), _matrix.rows()); i++)
    {
        singInvMatrix(i, i) = 1.0 / svd.singularValues()[i];
    }
    Eigen::MatrixXd U = svd.matrixU();
    Eigen::MatrixXd pinv = svd.matrixV()*singInvMatrix*U.transpose();
    return pinv;
}

static Eigen::MatrixXd weightedPinv(const Eigen::MatrixXd& _matrix, const Eigen::MatrixXd& _weight)
{
    Eigen::MatrixXd pinv;
    if (_matrix.rows() < _matrix.cols())
    {
        pinv = _weight.inverse() * _matrix.transpose() * (_matrix * _weight.inverse() * _matrix.transpose()).inverse();
    }
    else
    {
        pinv = (_matrix.transpose()*_matrix + _weight).inverse()*_matrix.transpose();
    }
    return pinv;
}

static Eigen::MatrixXd dampedLeastSquareInverse(const Eigen::MatrixXd& _matrix, const Eigen::MatrixXd& _weight)
{
    Eigen::MatrixXd pinv;
    if (_matrix.rows() < _matrix.cols())
    {
        Eigen::JacobiSVD<Eigen::MatrixXd> svd(_matrix, Eigen::ComputeFullU | Eigen::ComputeFullV);
        double lambda_sq;
        if (svd.singularValues()[_matrix.rows() - 1] > EPS_DLSI)
            lambda_sq = 0.0;
        else
        {
            double a = svd.singularValues()[_matrix.rows() - 1] / EPS_DLSI;
            lambda_sq = (1.0 - a*a)*lambdamax_DLSI*lambdamax_DLSI;
        }
        double beta = 1.0e-12 * lambda_sq;
        pinv = _weight.inverse() * _matrix.transpose() * (_matrix * _weight.inverse() * _matrix.transpose() + beta*Eigen::MatrixXd::Identity(_matrix.rows(), _matrix.rows())
                                                          + lambda_sq*svd.matrixU().col(_matrix.rows() - 1)*svd.matrixU().col(_matrix.rows() - 1).transpose()).inverse();
    }
    else
    {
        pinv = (_matrix.transpose()*_matrix + _weight).inverse()*_matrix.transpose();
    }
    return pinv;
}


static SO3 projectToSO3(Eigen::Matrix3d& R)
{
    double detR = R.deterstd::minant();
    Eigen::JacobiSVD<Eigen::Matrix3d> svd(R, Eigen::ComputeFullU | Eigen::ComputeFullV);

    if (detR > 0)
        R = svd.matrixU()*svd.matrixV().transpose();
    else
    {
        R = svd.matrixU().col(0)*svd.matrixV().col(0).transpose() + svd.matrixU().col(1)*svd.matrixV().col(1).transpose() - svd.matrixU().col(2)*svd.matrixV().col(2).transpose();
    }

    SO3 _R;
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
            _R(i, j) = R(i, j);
    }
    return _R;
}

static SE3 projectToSE3(SE3& T)
{
    Eigen::Matrix3d R;
    SO3 _R = T.GetOrientation();
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
            R(i, j) = _R(i, j);
    }
    T.SetOrientation(projectToSO3(R));
    return T;
}

static Eigen::VectorXd VecToVector(const std::vector<double>& vec)
{
    Eigen::VectorXd vector(vec.size());
    for (unsigned int i = 0; i < vec.size(); i++)
        vector(i) = vec[i];
    return vector;
}

static std::vector<double> VectorToVec(const Eigen::VectorXd& vect)
{
    std::vector<double> vec(vect.size());
    for (int i = 0; i < vect.size(); i++)
        vec[i] = vect(i);
    return vec;
}

static std::vector<double> SO3ToEulerZYX(const SO3& R)
{
    std::vector<double> euler(3);
    euler[1] = -asin(R[2]);
    euler[0] = atan2(R[1], R[0]);
    euler[2] = atan2(R[5], R[8]);
    return euler;
}

static std::vector<double> SO3ToEulerXYZ(const SO3& R)
{
    std::vector<double> euler(3);
    euler[1] = asin(R[6]);
    euler[0] = atan2(-R[7], R[8]);
    euler[2] = atan2(-R[3], R[0]);
    return euler;
}


static Vec3 VectortoVec3(const Eigen::Vector3d& _vector)
{
    Vec3 v(_vector(0), _vector(1), _vector(2));
    return v;
}

static Eigen::Matrix3d SO3toMatrix(const SO3& _SO3)
{
    Eigen::Matrix3d R;
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
            R(j, i) = _SO3(j, i);
    }
    return R;
}

static Eigen::Matrix4d SE3toMatrix(const SE3& _SE3)
{
    Eigen::Matrix4d T;
    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 3; j++)
            T(j, i) = _SE3(j, i);
    }
    T(3, 3) = 1;
    return T;
}

static Eigen::VectorXd SE3toVectorXd(const SE3& _SE3)
{
    Eigen::VectorXd vec(12);

    SO3 R = _SE3.GetOrientation();
    Vec3 p = _SE3.GetPosition();

    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            vec(i * 3 + j) = R(j, i);

    for (int i = 0; i < 3; i++)
        vec(i + 9) = p[i];

    return vec;

}

static SE3 VectorXdtoSE3(const Eigen::VectorXd vec)
{
    SO3 R;
    Vec3 p;

    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            R(j, i) = vec[i * 3 + j];

    for (int i = 0; i < 3; i++)
        p[i] = vec[i + 9];

    SE3 T(R, p);
    return T;
}

static Eigen::Vector3d Vec3toVector(const Vec3& _vec3)
{
    Eigen::Vector3d v(_vec3[0], _vec3[1], _vec3[2]);
    return v;
}

static Eigen::VectorXd se3toVector(const se3& _se3)
{
    Eigen::VectorXd v(6);
    for (int i = 0; i < 6; i++)
        v(i) = _se3[i];
    return v;
}

static Eigen::VectorXd dse3toVector(const dse3& _dse3)
{
    Eigen::VectorXd f(6);
    for (int i = 0; i < 6; i++)
        f(i) = _dse3[i];
    return f;
}

static dse3 Vectortodse3(const Eigen::VectorXd& _vector)
{
    dse3 _dse3;
    for (int i = 0; i < 6; i++)
        _dse3[i] = _vector(i);
    return _dse3;
}

static se3 Vectortose3(const Eigen::VectorXd& _vector)
{
    se3 _se3;
    for (int i = 0; i < 6; i++)
        _se3[i] = _vector(i);
    return _se3;
}

static Eigen::VectorXd Matrixse3toVectorse3(const Eigen::MatrixXd& _se3)
{
    Eigen::VectorXd __se3(6);
    __se3(0) = _se3(2, 1);
    __se3(1) = _se3(0, 2);
    __se3(2) = _se3(1, 0);
    for (int i = 0; i < 3; i++)
    {
        __se3(i + 3) = _se3(i, 3);
    }
    return __se3;
}

static Eigen::VectorXd concatenateVec3(const Vec3& vec1, const Vec3& vec2)
{
    Eigen::VectorXd vec(6);
    for (int i = 0; i < 3; i++) {
        vec(i) = vec1[i];
        vec(i + 3) = vec2[i];
    }
    return vec;
}

static Eigen::VectorXd concatenateVec6(const Eigen::VectorXd& vec1, const Eigen::VectorXd& vec2)
{
    Eigen::VectorXd vec(12);
    for (int i = 0; i < 6; i++) {
        vec(i) = vec1[i];
        vec(i +6) = vec2[i];
    }
    return vec;
}
static Eigen::VectorXd SE3toVector(const SE3& T)
{
    SO3 R = T.GetOrientation();
    Vec3 p = T.GetPosition();
    Vec3 r = Log(R);
    Eigen::VectorXd vec = concatenateVec3(r, p);
    return vec;
}

static SE3 VectortoSE3(const Eigen::VectorXd& vec)
{
    Vec3 r, p;
    for (int i = 0; i < 3; i++) {
        r[i] = vec(i);
        p[i] = vec(3 + i);
    }

    SO3 R = Exp(r);
    SE3 T(R, p);
    return T;
}

static Eigen::MatrixXd InertiatoMatrix(const Inertia& G)
{
    Eigen::MatrixXd M(6, 6);
    M(0, 0) = G[0];		M(0, 1) = G[3];		M(0, 2) = G[4];		M(0, 3) = 0.0;		M(0, 4) = -G[8];	M(0, 5) = G[7];
    M(1, 0) = G[3];		M(1, 1) = G[1];		M(1, 2) = G[5];		M(1, 3) = G[8];		M(1, 4) = 0.0;		M(1, 5) = -G[6];
    M(2, 0) = G[4];		M(2, 1) = G[5];		M(2, 2) = G[2];		M(2, 3) = -G[7];	M(2, 4) = G[6];		M(2, 5) = 0.0;
    M(3, 0) = 0.0;		M(3, 1) = G[8];		M(3, 2) = -G[7];	M(3, 3) = G[9];		M(3, 4) = 0.0;		M(3, 5) = 0.0;
    M(4, 0) = -G[8];	M(4, 1) = 0.0;		M(4, 2) = G[6];		M(4, 3) = 0.0;		M(4, 4) = G[9];		M(4, 5) = 0.0;
    M(5, 0) = G[7];		M(5, 1) = -G[6];	M(5, 2) = 0.0;		M(5, 3) = 0.0;		M(5, 4) = 0.0;		M(5, 5) = G[9];
    return M;
}

static Eigen::Matrix3d skewVec3(const Vec3& p)
{
    Eigen::Matrix3d skewP;
    skewP.setZero();
    skewP(0, 1) = -p[2];
    skewP(0, 2) = p[1];
    skewP(1, 0) = p[2];
    skewP(1, 2) = -p[0];
    skewP(2, 0) = -p[1];
    skewP(2, 1) = p[0];
    return skewP;
}

static Eigen::MatrixXd AdtoMatrix(const SE3& T)
{
    Eigen::MatrixXd AdjT(6, 6);
    AdjT.setZero();
    Eigen::Matrix3d Zero;
    Zero.setZero();
    AdjT << SO3toMatrix(T.GetOrientation()), Zero, skewVec3(T.GetPosition())*SO3toMatrix(T.GetOrientation()), SO3toMatrix(T.GetOrientation());
    return AdjT;
}

static Eigen::MatrixXd InvdAdtoMatrix(const SE3& T)
{
    Eigen::MatrixXd InvdAdjT(6, 6);
    InvdAdjT.setZero();
    Eigen::Matrix3d Zero;
    Zero.setZero();
    InvdAdjT << SO3toMatrix(T.GetOrientation()), skewVec3(T.GetPosition())*SO3toMatrix(T.GetOrientation()), Zero, SO3toMatrix(T.GetOrientation());
    return InvdAdjT;
}

static Axis Vec3toAxis(const Vec3& p)
{
    Axis a(p[0], p[1], p[2]);
    return a;
}

static SO3 MatrixtoSO3(const Eigen::Matrix3d& _R)
{
    // project _R to SO3
    Eigen::MatrixXd prjR;
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(_R, Eigen::ComputeFullU | Eigen::ComputeFullV);
    if (_R.deterstd::minant() > 0)
        prjR = svd.matrixU()*svd.matrixV().transpose();
    else {
        Eigen::Vector3d v(1, 1, -1);
        prjR = svd.matrixU()*v.asDiagonal()*svd.matrixV().transpose();
    }
    prjR(1, 1);
    SO3 R(prjR(0, 0), prjR(1, 0), prjR(2, 0), prjR(0, 1), prjR(1, 1), prjR(2, 1), prjR(0, 2), prjR(1, 2), prjR(2, 2));
    return R;
}

static SE3 SKKUtoSE3(const std::vector<double>& ori, const std::vector<double>& pos)
{
    return SE3(ori[0], ori[1], ori[2], ori[3], ori[4], ori[5], ori[6], ori[7], ori[8], pos[0], pos[1], pos[2]);
}

static double distSE3(const SE3 T1, const SE3 T2)
{
    double dist_p = Norm(T1.GetPosition() - T2.GetPosition());
    double dist_R = Norm(Log(Inv(T1.GetOrientation()) * T2.GetOrientation()));
    return dist_p + dist_R;
}

//
//static std::vector<pair<Vec3, SE3>> makeCylinderWithBoxes(SE3 cylinderCenter, double cylinderRadius, double cylinderHeight, double thickness, double epsilon, int numBox)
//{
//	// axis of along cylinder height must be aligned to z-axis
//	// numBox must be larger than 2
//
//	std::vector<pair<Vec3, SE3>> boxGeostd::minfo(numBox);
//
//	double tmpR = cylinderRadius - thickness / 2;
//	double tmpTheta;
//	Vec3 tmpVec3;
//	SE3 tmpSE3;
//	double axisVal1, axisVal2;
//
//
//	tmpTheta = 2 * SR_PI * 0.5 / ((double)numBox);
//
//	// if numBox is larger than 2, axisVal1 and axisVal2 have positive values
//	axisVal1 = sqrt((tmpR*tmpR) / (1 + tan(tmpTheta)*tan(tmpTheta)));
//	axisVal2 = axisVal1 * tan(tmpTheta);
//
//	tmpVec3[0] = thickness;
//	tmpVec3[1] = 2 * ((cylinderRadius - thickness) * tan(tmpTheta) - epsilon);
//	tmpVec3[2] = cylinderHeight;
//	tmpSE3 = RotZ(tmpTheta);
//	tmpSE3.SetPosition(Vec3(tmpR*cos(tmpTheta), tmpR*sin(tmpTheta), 0.0));
//	tmpSE3 = tmpSE3;
//
//	boxGeostd::minfo[0].first = tmpVec3;
//	boxGeostd::minfo[0].second = tmpSE3;
//
//	for (int i = 1; i < numBox; i++)
//	{
//		boxGeostd::minfo[i].first = tmpVec3;
//		tmpSE3 = RotZ(2 * SR_PI / ((double)numBox)) * tmpSE3;
//		boxGeostd::minfo[i].second = tmpSE3;
//	}
//
//	for (int i = 0; i < numBox; i++)
//		boxGeostd::minfo[i].second = cylinderCenter * boxGeostd::minfo[i].second;
//
//
//	return boxGeostd::minfo;
//}
//
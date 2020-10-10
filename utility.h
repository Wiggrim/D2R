#pragma once

#include <math.h>
#include <eigen3/Eigen/Dense>

#define EM_a			6378137.0
#define EM_b			6356752.3142
#define EM_e2			(1-(EM_b*EM_b/EM_a/EM_a))
#define PI				3.1415926535897932
#define PAR_C			2.99792458e8
#define EM_we           7.292e-5

namespace Utility
{
    double CalcGravity(double latitude, double hae);
    bool MatrixMultipy(double *c, double *a, int dim_am, int dim_an, double *b, int dim_bm, int dim_bn);
    void CalcR(double *M, double *N, double a, double e2, double lat);
    void Pos2ECEF(double *ecef_xyz, double e2, double Rn, double lat, double lon, double hae);
    void Pos2DCM(double *C_ne, double lat, double lon);
    bool Dcm2Euler( Eigen::Matrix3d& C, Eigen::Vector3d& Euler );
    void SkewMatrix( Eigen::Vector3d& vec, Eigen::Matrix3d& skew_mat );
}
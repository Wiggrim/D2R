#pragma once

#include <eigen3/Eigen/Dense>
#include "parameters.h"

using namespace Eigen;

class IntegrationBase
{
  public:
    IntegrationBase() = delete;
    IntegrationBase(const Eigen::Vector3d &_acc_0, const Eigen::Vector3d &_gyr_0,
                    const Eigen::Vector3d &_linearized_ba, const Eigen::Vector3d &_linearized_bg)
        : acc_0{_acc_0}, gyr_0{_gyr_0}, linearized_acc{_acc_0}, linearized_gyr{_gyr_0},
          linearized_ba{_linearized_ba}, linearized_bg{_linearized_bg},
            jacobian{Eigen::Matrix<double, 15, 15>::Identity()}, covariance{Eigen::Matrix<double, 15, 15>::Zero()},
          sum_dt{0.0}, delta_p{Eigen::Vector3d::Zero()}, delta_q{Eigen::Quaterniond::Identity()}, delta_v{Eigen::Vector3d::Zero()}

    {
        noise = Eigen::Matrix<double, 12, 12>::Zero();
        noise.block( 0, 0, 3, 3 ) = Eigen::Matrix3d::Identity() * ACC_NOISE_STD * ACC_NOISE_STD;
        noise.block( 3, 3, 3, 3 ) = Eigen::Matrix3d::Identity() * GYR_NOISE_STD * GYR_NOISE_STD;
        noise.block( 6, 6, 3, 3 ) = Eigen::Matrix3d::Identity() * ACC_DRIFT_STD * ACC_DRIFT_STD;
        noise.block( 9, 9, 3, 3 ) = Eigen::Matrix3d::Identity() * GYR_DRIFT_STD * GYR_DRIFT_STD;
    }

    void push_back(double dt, const Eigen::Vector3d &acc, const Eigen::Vector3d &gyr)
    {
        dt_buf.push_back(dt);
        acc_buf.push_back(acc);
        gyr_buf.push_back(gyr);
        propagate2(dt, acc, gyr);
    }

    void midPointIntegration2( double _dt_0, double _dt_1,
                               const Eigen::Vector3d &_acc_0, const Eigen::Vector3d &_gyr_0,
                               const Eigen::Vector3d &_acc_1, const Eigen::Vector3d &_gyr_1,
                               const Eigen::Vector3d &delta_p, const Eigen::Quaterniond &delta_q, const Eigen::Vector3d &delta_v,
                               const Eigen::Vector3d &linearized_ba, const Eigen::Vector3d &linearized_bg,
                               Eigen::Vector3d &result_delta_p, Eigen::Quaterniond &result_delta_q, Eigen::Vector3d &result_delta_v,
                               Eigen::Vector3d &result_linearized_ba, Eigen::Vector3d &result_linearized_bg, bool update_jacobian
                             )
    {
        Vector3d un_acc_0 = _dt_0 * ( _acc_0 - linearized_ba );
        Vector3d un_acc_1 = _dt_1 * ( _acc_1 - linearized_ba );
        Vector3d un_gyr_0 = _dt_0 * ( _gyr_0 - linearized_bg );
        Vector3d un_gyr_1 = _dt_1 * ( _gyr_1 - linearized_bg );
        
        Vector3d acc_rot = 0.5 * un_gyr_1.cross( un_acc_1 );
        Vector3d acc_scul = ( un_gyr_0.cross( un_acc_1 ) - un_gyr_1.cross( un_acc_0 ) ) / 12.0;
        Vector3d gyr_con = ( un_gyr_0.cross( un_gyr_1 ) ) / 12.0;
        
        Vector3d un_acc = un_acc_1 + acc_rot + acc_scul;
        Vector3d un_gyr = un_gyr_1 + gyr_con;
        
        result_delta_v = delta_v + delta_q * un_acc;
        result_delta_q = delta_q * Quaterniond( 1, un_gyr(0)/2, un_gyr(1)/2, un_gyr(2)/2 );
        result_delta_p = delta_p + ( delta_v + result_delta_v ) / 2.0 * _dt_1;
        
        result_linearized_ba = linearized_ba;
        result_linearized_bg = linearized_bg;   
        
        if( update_jacobian )
        {
            Matrix3d cross_gyr;
            cross_gyr << 0, -un_gyr_1(2), un_gyr_1(1), 
                        un_gyr_1(2), 0, -un_gyr_1(0),
                        -un_gyr_1(1), un_gyr_1(0), 0;  
                        
            Matrix3d cross_acc;
            cross_acc << 0, -un_acc_1(2), un_acc_1(1),
                        un_acc_1(2), 0, -un_acc_1(0),
                        -un_acc_1(1), un_acc_1(0), 0;
            
            MatrixXd F = MatrixXd::Zero( 15, 15 );
            
            F.block<3,3>( 0,0 ) = Matrix3d::Identity();
            F.block<3,3>( 0,6 ) = Matrix3d::Identity() * _dt_1;
            
            F.block<3,3>( 3,3 ) = Matrix3d::Identity() - cross_gyr;
            F.block<3,3>( 3,12 ) = -Matrix3d::Identity() * _dt_1;
            
            F.block<3,3>( 6,3 ) = -delta_q.toRotationMatrix() * cross_acc;
	          F.block<3,3>( 6,6 ) = Matrix3d::Identity();
            F.block<3,3>( 6,9 ) = -delta_q.toRotationMatrix() * _dt_1;
            
            F.block<3,3>( 9,9 ) = Matrix3d::Identity();
            
            F.block<3,3>( 12,12 ) = Matrix3d::Identity();
            
            MatrixXd V = MatrixXd::Zero( 15, 12 );
            
            V.block<3,3>( 3, 3 ) = -Matrix3d::Identity();
            
            V.block<3,3>( 6, 0 ) = -delta_q.toRotationMatrix();
            
            V.block<3,3>( 9, 6 ) = Matrix3d::Identity();
            
            V.block<3,3>( 12, 9 ) = Matrix3d::Identity();
            
            jacobian = F * jacobian;
            covariance = F * covariance * F.transpose() + V * noise * V.transpose() * _dt_1;            
        }
    }

    void propagate2(double _dt, const Eigen::Vector3d &_acc_1, const Eigen::Vector3d &_gyr_1)
    {
        double dt0;
        dt0 = dt;
        dt = _dt;
        acc_1 = _acc_1;
        gyr_1 = _gyr_1;
        Vector3d result_delta_p;
        Quaterniond result_delta_q;
        Vector3d result_delta_v;
        Vector3d result_linearized_ba;
        Vector3d result_linearized_bg;

        midPointIntegration2(dt0, _dt, acc_0, gyr_0, _acc_1, _gyr_1, delta_p, delta_q, delta_v,
                            linearized_ba, linearized_bg,
                            result_delta_p, result_delta_q, result_delta_v,
                            result_linearized_ba, result_linearized_bg, 1);

        delta_p = result_delta_p;
        delta_q = result_delta_q;
        delta_v = result_delta_v;
        linearized_ba = result_linearized_ba;
        linearized_bg = result_linearized_bg;
        delta_q.normalize();
        sum_dt += dt;
        acc_0 = acc_1;
        gyr_0 = gyr_1;  
     
    }

    double dt;
    Eigen::Vector3d acc_0, gyr_0;
    Eigen::Vector3d acc_1, gyr_1;

    const Eigen::Vector3d linearized_acc, linearized_gyr;
    Eigen::Vector3d linearized_ba, linearized_bg;

    Eigen::Matrix<double, 15, 15> jacobian, covariance;
    Eigen::Matrix<double, 15, 15> step_jacobian;
    Eigen::Matrix<double, 15, 18> step_V;
    Eigen::Matrix<double, 12, 12> noise;

    double sum_dt;
    Eigen::Vector3d delta_p;
    Eigen::Quaterniond delta_q;
    Eigen::Vector3d delta_v;

    std::vector<double> dt_buf;
    std::vector<Eigen::Vector3d> acc_buf;
    std::vector<Eigen::Vector3d> gyr_buf;

};

#include "Dop2Rot.h"

namespace D2R
{
    Filter::Filter()
    {
        initial_flag = false;

        cur_rot = Eigen::Quaterniond::Identity();
        cur_gyr_bias = Eigen::Vector3d::Zero();
        covariance = Eigen::MatrixXd::Zero( 6, 6 );
        acc_pre = Eigen::Vector3d::Zero();
        gyr_pre = Eigen::Vector3d::Zero();
        dt_pre = 0.;
        run_time = 0.;
    }

    bool Filter::Initialize( Eigen::Quaterniond& init_rot, Eigen::Matrix3d& init_rot_cov, Eigen::Vector3d& acc_pre, Eigen::Vector3d& gyr_pre )
    {
        cur_rot = init_rot;
        cur_gyr_bias = Eigen::Vector3d::Zero();
        covariance.block<3,3>( 0, 0 ) = init_rot_cov;
        covariance.block<3,3>( 3, 3 ) = Eigen::Matrix3d::Identity() * INIT_GYR_STD * INIT_GYR_STD;
        this->acc_pre = acc_pre;
        this->gyr_pre = gyr_pre;
        dt_pre = 0.01;

        initial_flag = true;

        return true;
    }

    bool Filter::PushBack( IMUBatches& imu_batches )
    {
        if( initial_flag == false )
        {
            std::cout << "Push Back called without initialization\n";
            return false;
        }

        unsigned int imu_meas_num = imu_batches.size();
        if( imu_meas_num <= 0 ) 
        {
            std::cout << "Too little IMU number for Filtering\n";
            return false;
        }

        for( IMUMeas& imu_meas : imu_batches )
        {
            Eigen::Vector3d acc_cur = imu_meas.acc;
            Eigen::Vector3d gyr_cur = imu_meas.gyr;
            double dt_cur = imu_meas.dt;

            Eigen::Vector3d un_rot_0 = ( gyr_pre - cur_gyr_bias ) * dt_pre;
            Eigen::Vector3d un_rot_1 = ( gyr_cur - cur_gyr_bias ) * dt_cur;
            Eigen::Vector3d rot_con = un_rot_0.cross( un_rot_1 ) / 12.;

            Eigen::Vector3d un_rot = un_rot_1 + rot_con;
            Eigen::Matrix3d skew_rot;
            skew_rot << 0., -un_rot(2), un_rot(1), \
                                    un_rot(2), 0., -un_rot(0), \
                                    -un_rot(1), un_rot(0), 0.;

            cur_rot = cur_rot * Eigen::Quaterniond( 1., un_rot.x() / 2., un_rot.y() / 2., un_rot.z() / 2. );
            cur_rot.normalize();

            cur_gyr_bias = cur_gyr_bias;

            Eigen::Matrix<double, 6, 6> F = Eigen::MatrixXd::Zero( 6, 6 );
            F.block<3,3>( 0, 0 ) = Eigen::Matrix3d::Identity() - skew_rot;
            F.block<3,3>( 0, 3 ) = -Eigen::Matrix3d::Identity() * dt_cur;
            F.block<3,3>( 3, 3 ) = Eigen::Matrix3d::Identity();

            Eigen::Matrix<double, 6, 6> Q = Eigen::MatrixXd::Zero( 6, 6 );
            Q.block<3,3>( 0, 0 ) = Eigen::Matrix3d::Identity() * GYR_NOISE_STD * GYR_NOISE_STD * dt_cur;
            Q.block<3,3>( 3, 3 ) = Eigen::Matrix3d::Identity() * GYR_DRIFT_STD * GYR_DRIFT_STD * dt_cur;

            covariance = F * covariance * F.transpose() + Q;

            acc_pre = acc_cur;
            gyr_pre = gyr_cur;
            dt_pre = dt_cur;
        }

        acc_pre = imu_batches[imu_meas_num-1].acc;
        gyr_pre = imu_batches[imu_meas_num-1].gyr;
        dt_pre = imu_batches[imu_meas_num-1].dt;

        return true;
    }

    bool Filter::Update( std::vector<DeltaDop>& delta_dops, IntegrationBase* integration_ptr )
    {
        if( initial_flag == false )
        {
            std::cout << "Filter called without initialization\n";
            return false;
        }

        if( integration_ptr == nullptr )
        {
            std::cout << "Invalid Integration Base Pointer\n";
            return false;
        }

        if( integration_ptr->sum_dt == 0. )
        {
            std::cout << "Invalid Integration Base Data\n";
            return false;
        }

        Eigen::Matrix3d gyr_2_delta_v = integration_ptr->jacobian.block( 6, 12 ,3, 3 );
        Eigen::Vector3d calibd_delta_v = integration_ptr->delta_v + gyr_2_delta_v * cur_gyr_bias;
        Eigen::Matrix3d skew_calibd_delta_v;
        Utility::SkewMatrix( calibd_delta_v, skew_calibd_delta_v );

        unsigned int sat_meas_num = delta_dops.size();
        if( sat_meas_num <= 0 )
        {
            std::cout << "There is no data for update\n";
            return true;
        }

        Eigen::MatrixXd los_mat = Eigen::MatrixXd::Zero( sat_meas_num, 3 );
        Eigen::MatrixXd J = Eigen::MatrixXd::Zero( sat_meas_num, 6 );
        Eigen::MatrixXd Inno = Eigen::MatrixXd::Zero( sat_meas_num, 1 );
        Eigen::MatrixXd Q = Eigen::MatrixXd::Ones( sat_meas_num, sat_meas_num ) * OSC_STD * OSC_STD + Eigen::MatrixXd::Identity( sat_meas_num, sat_meas_num ) * DOP_STD * DOP_STD;
        Eigen::MatrixXd K = Eigen::MatrixXd::Zero( 6, sat_meas_num );
        Eigen::Matrix<double, 6, 6> P = covariance;
        for( unsigned int i = 0; i < sat_meas_num; i++ )
        {
            Eigen::Vector3d los = delta_dops[i].los; 
            double delta_d = delta_dops[i].deltaD; 
            los_mat.block( i, 0, 1, 3 ) = los.transpose();

            Inno( i, 0 ) = delta_d - los.transpose() * cur_rot.toRotationMatrix() * calibd_delta_v;
            J.block( i, 0, 1, 3 ) = -los.transpose() * cur_rot.toRotationMatrix() * skew_calibd_delta_v;
            J.block( i, 3, 1, 3 ) = los.transpose() * cur_rot.toRotationMatrix() * gyr_2_delta_v;
        }

        Eigen::MatrixXd covariance_obs_inv =  ( J*P*J.transpose() + Q ).inverse();
        K = P * J.transpose() * covariance_obs_inv;

        Eigen::Matrix<double, 6, 1> state_update = K * Inno;
        cur_rot = cur_rot * Eigen::Quaterniond( 1., state_update(0,0) / 2., state_update(1,0) / 2., state_update(2,0) / 2. );
        cur_rot.normalize();
        cur_gyr_bias = cur_gyr_bias + state_update.block( 3, 0, 3, 1 );        

        P = P - K * J * P;
        P = (P + P.transpose()) / 2.;
        covariance = P;

        run_time += integration_ptr->sum_dt;

        return true;
    }


}
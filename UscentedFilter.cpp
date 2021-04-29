#include "Dop2Rot.h"

namespace D2R
{
    UnscentedFilter::UnscentedFilter()
    {
        initial_flag = false;
        cur_rot.setIdentity();
        cur_gyr_bias.setZero();
        covariance.setZero();
        run_time = 0.;
        acc_pre.setZero();
        gyr_pre.setZero();
        dt_pre = 0.;
    }

    bool UnscentedFilter::Initialize( Eigen::Quaterniond& init_rot, Eigen::Matrix3d& init_rot_cov, Eigen::Vector3d& acc_pre, Eigen::Vector3d& gyr_pre )
    {
        cur_rot = init_rot;
        cur_gyr_bias.setZero();
        covariance.block<3,3>( 0, 0 ) = init_rot_cov;
        covariance.block<3,3>( 3, 3 ) = Eigen::Matrix3d::Identity() * INIT_GYR_STD * INIT_GYR_STD;
        this->acc_pre = acc_pre;
        this->gyr_pre = gyr_pre;
        dt_pre = 0.01;
        initial_flag = true;

        return true;
    }

    bool UnscentedFilter::PushBack( IMUBatches& imu_batches )
    {
        if( initial_flag == false )
        {
            std::cout << "Should called after initialized\n";
            return false;
        }

        if( imu_batches.size() <= 0 )
        {
            std::cout << "Little values in IMU measurements batches\n";
            return false;
        }

        Eigen::Matrix<double, 6, 6> L = covariance.llt().matrixL();
        double k = 3;
        double n = 6;
        double scale = sqrt( n+k );
        L = L * scale;

        std::vector<Eigen::Quaterniond, Eigen::aligned_allocator<Eigen::Quaterniond>> input_sigma_points_rot;
        std::vector<Eigen::Vector3d, Eigen::aligned_allocator<Eigen::Vector3d>> input_sigma_points_bias;
        input_sigma_points_rot.resize( 13 );
        input_sigma_points_bias.resize( 13 );
        for( unsigned int i = 0; i < 6; i++ )
        {
            input_sigma_points_rot[i] = cur_rot * Eigen::Quaterniond( 1., L(0,i)/2., L(1,i)/2., L(2,i)/2. );
            input_sigma_points_rot[i].normalize();
            input_sigma_points_bias[i] = cur_gyr_bias + L.block( 3, i, 3, 1 );
        }
        for( unsigned int i = 0; i < 6; i++ )
        {
            input_sigma_points_rot[i+6] = cur_rot * Eigen::Quaterniond( 1., -L(0,i)/2., -L(1,i)/2., -L(2,i)/2. );
            input_sigma_points_rot[i+6].normalize();
            input_sigma_points_bias[i+6] = cur_gyr_bias - L.block( 3, i, 3, 1 );            
        }
        input_sigma_points_rot[12] = cur_rot;
        input_sigma_points_bias[12] = cur_gyr_bias;

        Eigen::Matrix4d cov_q = Eigen::Matrix4d::Zero();
        Eigen::Vector3d acc_pre_latch;
        Eigen::Vector3d gyr_pre_latch;
        double dt_pre_latch;
        for( unsigned int i = 0; i < 13; i++ )
        {
            acc_pre_latch = acc_pre;
            gyr_pre_latch = gyr_pre;
            dt_pre_latch = dt_pre;
            Eigen::Vector3d cur_gyr_bias_latch = input_sigma_points_bias[i];
            Eigen::Quaterniond cur_rot_latch = input_sigma_points_rot[i];
            for( IMUMeas& imu_meas : imu_batches )
            {
                Eigen::Vector3d acc_cur = imu_meas.acc;
                Eigen::Vector3d gyr_cur = imu_meas.gyr;
                double dt_cur = imu_meas.dt;

                Eigen::Vector3d un_rot_0 = ( gyr_pre_latch - cur_gyr_bias_latch ) * dt_pre_latch;
                Eigen::Vector3d un_rot_1 = ( gyr_cur - cur_gyr_bias_latch ) * dt_cur;
                Eigen::Vector3d rot_con = un_rot_0.cross( un_rot_1 ) / 12.;

                Eigen::Vector3d un_rot = un_rot_1 + rot_con;
                cur_rot_latch = cur_rot_latch * Eigen::Quaterniond( 1., un_rot.x() / 2., un_rot.y() / 2., un_rot.z() / 2. );
                cur_rot_latch.normalize();

                acc_pre_latch = acc_cur;
                gyr_pre_latch = gyr_cur;
                dt_pre_latch = dt_cur;
            }
            input_sigma_points_rot[i] = cur_rot_latch;
            input_sigma_points_bias[i] = cur_gyr_bias_latch;
            Eigen::Vector4d vec_q = Eigen::Vector4d( cur_rot_latch.x(), cur_rot_latch.y(), cur_rot_latch.z(), cur_rot_latch.w() );

            if( i != 12 )
            {
                cov_q = cov_q + vec_q * vec_q.transpose() / 2. / (n+k);
            }
            else
            {
                cov_q = cov_q + vec_q * vec_q.transpose() / 3.;
            }
            
        }

        Eigen::SelfAdjointEigenSolver<Eigen::Matrix4d> eigen_cov_q( cov_q, Eigen::ComputeEigenvectors );
        Eigen::Quaterniond mean_rot = Eigen::Quaterniond( eigen_cov_q.eigenvectors()(3,3), eigen_cov_q.eigenvectors()(0,3), eigen_cov_q.eigenvectors()(1,3), eigen_cov_q.eigenvectors()(2,3) );
        Eigen::Vector3d mean_bias = cur_gyr_bias;
        Eigen::Matrix<double, 6, 6> covariance_next = Eigen::MatrixXd::Zero( 6, 6 );
        for( unsigned int i = 0; i < 13; i++ )
        {
            Eigen::Matrix<double, 6, 1> delta_vec;
            delta_vec.setZero();

            Eigen::Quaterniond delta_rot = mean_rot.inverse() * input_sigma_points_rot[i];
            Eigen::Vector3d delta_bias = input_sigma_points_bias[i] - mean_bias;
            Eigen::Vector3d delta_rot_vec;
            if( delta_rot.w() >= 0. ) delta_rot_vec = delta_rot.vec() * 2.;
            else delta_rot_vec = delta_rot.vec() * -2.;

            delta_vec.block<3,1>( 0, 0 ) = delta_rot_vec;
            delta_vec.block<3,1>( 3, 0 ) = delta_bias;

            if( i != 12 )
            {
                covariance_next += delta_vec * delta_vec.transpose() / 2. / (n+k);
            }
            else
            {
                covariance_next += delta_vec * delta_vec.transpose() / 3.;
            }
        }

        double sum_t = 0.;
        for( unsigned int i = 0; i < imu_batches.size(); i++ )
        {
            sum_t += imu_batches[i].dt;
        }

        Eigen::Matrix<double, 6, 6> Q = Eigen::MatrixXd::Zero( 6, 6 );
        Q.block( 0, 0, 3, 3 ) = Eigen::Matrix3d::Identity() * GYR_NOISE_STD * GYR_NOISE_STD * sum_t;
        Q.block( 3, 3, 3, 3 ) = Eigen::Matrix3d::Identity() * GYR_DRIFT_STD * GYR_DRIFT_STD * sum_t;

        cur_rot = mean_rot;
        cur_gyr_bias = cur_gyr_bias;
        covariance = covariance_next + Q;
        acc_pre = imu_batches.back().acc;
        gyr_pre = imu_batches.back().gyr;
        dt_pre = imu_batches.back().dt;
        run_time += sum_t;

        return true;
    }

    bool UnscentedFilter::Update( Eigen::Vector3d& delta_v_local, IntegrationBase* integration_ptr, Eigen::Matrix3d& delta_v_local_cov )
    {
        if( integration_ptr == nullptr )
        {
            std::cout << "Invalid IMU integration pointer\n";
            return false;
        }

        if( initial_flag == false )
        {
            std::cout << "Update should be called after initialized\n";
            return false;
        }

        Eigen::Matrix<double, 6, 6> L = covariance.llt().matrixL();
        double n = 6;
        double k = 3;
        double scale = sqrt( n + k );
        L = L * scale;

        std::vector<Eigen::Quaterniond, Eigen::aligned_allocator<Eigen::Quaterniond>> input_sigma_points_rot;
        std::vector<Eigen::Vector3d, Eigen::aligned_allocator<Eigen::Vector3d>> input_sigma_points_bias;
        input_sigma_points_rot.resize( 13 );
        input_sigma_points_bias.resize( 13 );
        for( unsigned int i = 0; i < 6; i++ )
        {
            input_sigma_points_rot[i] = cur_rot * Eigen::Quaterniond( 1., L(0,i)/2., L(1,i)/2., L(2,i)/2. );
            input_sigma_points_rot[i].normalize();
            input_sigma_points_bias[i] = cur_gyr_bias + L.block( 3, i, 3, 1 );
        }
        for( unsigned int i = 0; i < 6; i++ )
        {
            input_sigma_points_rot[i+6] = cur_rot * Eigen::Quaterniond( 1., -L(0,i)/2., -L(1,i)/2., -L(2,i)/2. );
            input_sigma_points_rot[i+6].normalize();
            input_sigma_points_bias[i+6] = cur_gyr_bias - L.block( 3, i, 3, 1 );            
        }
        input_sigma_points_rot[12] = cur_rot;
        input_sigma_points_bias[12] = cur_gyr_bias;        

        Eigen::Matrix3d gyr_delta_v = integration_ptr->jacobian.block( 6, 12, 3, 3 );
        Eigen::Vector3d delta_v_b = integration_ptr->delta_v;
        std::vector<Eigen::Vector3d, Eigen::aligned_allocator<Eigen::Vector3d>> output_sigma_points_delta_v;
        Eigen::Vector3d mean_delta_v_local = Eigen::Vector3d::Zero();
        output_sigma_points_delta_v.resize( 13 );
        for( unsigned int i = 0; i < 13; i++ )
        {
            output_sigma_points_delta_v[i] = input_sigma_points_rot[i].toRotationMatrix() * (delta_v_b + gyr_delta_v * (input_sigma_points_bias[i] - integration_ptr->linearized_bg));
            if( i == 12 )
            {
                mean_delta_v_local += output_sigma_points_delta_v[i] / 3.;
            }
            else
            {
                mean_delta_v_local += output_sigma_points_delta_v[i] / 2. / (n+k);
            }
        }

        Eigen::Matrix3d cov_obs = Eigen::Matrix3d::Zero();
        Eigen::Matrix<double, 6, 3> cross_cov = Eigen::Matrix<double, 6, 3>::Zero();
        for( unsigned int i = 0; i < 6; i++ )
        {
            cov_obs += ( output_sigma_points_delta_v[i] - mean_delta_v_local ) * ( output_sigma_points_delta_v[i] - mean_delta_v_local ).transpose() / 2. / (n+k);
            cross_cov += L.block( 0, i, 6, 1 ) * ( output_sigma_points_delta_v[i] - mean_delta_v_local ).transpose() / 2. / (n+k);
        }
        for( unsigned int i = 0; i < 6; i++ )
        {
            cov_obs += ( output_sigma_points_delta_v[i+6] - mean_delta_v_local ) * ( output_sigma_points_delta_v[i+6] - mean_delta_v_local ).transpose() / 2. / (n+k);
            cross_cov += -L.block( 0, i, 6, 1 ) * ( output_sigma_points_delta_v[i+6] - mean_delta_v_local ).transpose() / 2. / (n+k);
        }
        cov_obs += ( output_sigma_points_delta_v[12] - mean_delta_v_local ) * ( output_sigma_points_delta_v[12] - mean_delta_v_local ).transpose() / 3.;
        cov_obs += delta_v_local_cov;

        Eigen::Vector3d delta_delta_v_local = delta_v_local - mean_delta_v_local;

        Eigen::Matrix<double, 6, 1> delta_state = cross_cov * cov_obs.inverse() * delta_delta_v_local;
        cur_rot = cur_rot * Eigen::Quaterniond( 1., delta_state(0) / 2., delta_state(1) / 2., delta_state(2) / 2. );
        cur_rot.normalize();
        cur_gyr_bias = cur_gyr_bias + delta_state.block( 3, 0, 3, 1 );
        covariance = covariance - cross_cov * cov_obs.inverse() * cross_cov.transpose();
        covariance = ( covariance + covariance.transpose() ) / 2.;

        return true;

    }

    bool UnscentedFilter::Update( std::vector<DeltaDop>& delta_dops, IntegrationBase* integration_ptr )
    {
        if( initial_flag == false )
        {
            std::cout << "UKF Update should be called after initialize\n";
            return false;
        }

        if( integration_ptr == nullptr )
        {
            std::cout << "IMU Integration in UKF Update should not be NULL\n";
            return false;
        }

        if( delta_dops.size() <= 0 )
        {
            std::cout << "Little GNSS Measurement in UKF Update\n";
            return false;
        }

        unsigned int sat_meas_num = delta_dops.size();
        Eigen::MatrixXd delta_d = Eigen::MatrixXd::Zero( sat_meas_num, 1 );
        Eigen::MatrixXd los_mat = Eigen::MatrixXd::Zero( sat_meas_num, 3 );
        for( unsigned int i = 0; i < sat_meas_num; i++ )
        {
            delta_d( i, 0 ) = delta_dops[i].deltaD;
            los_mat.block( i, 0, 1, 3 ) = delta_dops[i].los.transpose();
        }

        Eigen::Matrix<double, 6, 6> L = covariance.llt().matrixL();
        double n = 6;
        double k = 3;
        double scale = sqrt( n + k );
        L = L * scale;

        std::vector<Eigen::Quaterniond, Eigen::aligned_allocator<Eigen::Quaterniond>> input_sigma_points_rot;
        std::vector<Eigen::Vector3d, Eigen::aligned_allocator<Eigen::Vector3d>> input_sigma_points_bias;
        input_sigma_points_rot.resize( 13 );
        input_sigma_points_bias.resize( 13 );
        for( unsigned int i = 0; i < 6; i++ )
        {
            input_sigma_points_rot[i] = cur_rot * Eigen::Quaterniond( 1., L(0,i)/2., L(1,i)/2., L(2,i)/2. );
            input_sigma_points_rot[i].normalize();
            input_sigma_points_bias[i] = cur_gyr_bias + L.block( 3, i, 3, 1 );
        }
        for( unsigned int i = 0; i < 6; i++ )
        {
            input_sigma_points_rot[i+6] = cur_rot * Eigen::Quaterniond( 1., -L(0,i)/2., -L(1,i)/2., -L(2,i)/2. );
            input_sigma_points_rot[i+6].normalize();
            input_sigma_points_bias[i+6] = cur_gyr_bias - L.block( 3, i, 3, 1 );            
        }
        input_sigma_points_rot[12] = cur_rot;
        input_sigma_points_bias[12] = cur_gyr_bias;                

        Eigen::Matrix3d gyr_delta_v = integration_ptr->jacobian.block( 6, 12, 3, 3 );
        Eigen::Vector3d delta_v_b = integration_ptr->delta_v;
        Eigen::MatrixXd output_sigma_points_delta_d = Eigen::MatrixXd::Zero( sat_meas_num, 13 );
        Eigen::MatrixXd delta_d_mean = Eigen::MatrixXd::Zero( sat_meas_num, 1 );
        Eigen::MatrixXd delta_d_cov = Eigen::MatrixXd::Zero( sat_meas_num, sat_meas_num );
        Eigen::MatrixXd cross_cov = Eigen::MatrixXd::Zero( 6, sat_meas_num );
        for( unsigned int i = 0; i < 13; i++ )
        {
            output_sigma_points_delta_d.block( 0, i, sat_meas_num, 1 ) = los_mat * input_sigma_points_rot[i].toRotationMatrix() * \
                                                                                                                                                ( delta_v_b + gyr_delta_v * (input_sigma_points_bias[i] - integration_ptr->linearized_bg) );
            if( i != 12 ) delta_d_mean += output_sigma_points_delta_d.block( 0, i, sat_meas_num, 1 ) / 2. / (n+k);
            else delta_d_mean += output_sigma_points_delta_d.block( 0, i, sat_meas_num, 1 ) / 3.;
        }

        for( unsigned int i = 0; i < 6; i++ )
        {
            Eigen::MatrixXd delta = output_sigma_points_delta_d.block( 0, i, sat_meas_num, 1 ) - delta_d_mean;
            delta_d_cov += delta * delta.transpose() / 2. / (n+k);
            cross_cov += L.block( 0, i, 6, 1 ) * delta.transpose() / 2. / (n+k);
        }
        for( unsigned int i = 0; i < 6; i++ )
        {
            Eigen::MatrixXd delta = output_sigma_points_delta_d.block( 0, i+6, sat_meas_num, 1 ) - delta_d_mean;
            delta_d_cov += delta * delta.transpose() / 2. / (n+k);
            cross_cov += -L.block( 0, i, 6, 1 ) * delta.transpose() / 2. / (n+k);
        }
        {
            Eigen::MatrixXd delta = output_sigma_points_delta_d.block( 0, 12, sat_meas_num, 1 ) - delta_d_mean;
            delta_d_cov += delta * delta.transpose() / 3.;
        }
        Eigen::MatrixXd R = Eigen::MatrixXd::Identity( sat_meas_num, sat_meas_num ) * DOP_STD * DOP_STD * 2 + Eigen::MatrixXd::Ones( sat_meas_num, sat_meas_num ) * OSC_STD * OSC_STD;
        delta_d_cov += R;

        Eigen::MatrixXd Inno = delta_d - delta_d_mean;

        Eigen::MatrixXd K = cross_cov * delta_d_cov.inverse();
        Eigen::Matrix<double, 6, 1> state_update = K * Inno;
        cur_rot = cur_rot * Eigen::Quaterniond( 1., state_update(0) / 2., state_update(1) / 2., state_update(2) / 2. );
        cur_rot.normalize();
        cur_gyr_bias = cur_gyr_bias + state_update.block( 3, 0, 3, 1 );
        covariance = covariance - K * cross_cov.transpose();
        covariance = ( covariance + covariance.transpose() ) / 2.;

        return true;
    }

}
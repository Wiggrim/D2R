#pragma once

#include <eigen3/Eigen/Dense>
#include <vector>
#include <queue>
#include <algorithm>
#include <iostream>
#include <thread>
#include "integration_base.h"
#include "utility.h"
#include "parameters.h"

namespace D2R
{
    class SatData
    {
        public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW 

        unsigned int satNum;
        double transmitTime;
        double carrierDoppler;
        double carrierFrequency;
        Eigen::Vector3d los;
        Eigen::Vector3d satPos;
        Eigen::Vector3d satVel;

        bool operator== ( const SatData& sat_data ) { return this->satNum == sat_data.satNum; }
    };

    class IMUMeas
    {
        public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW 

        Eigen::Vector3d acc;
        Eigen::Vector3d gyr;
        double dt;
    };

    typedef std::vector<IMUMeas> IMUBatches;
    typedef std::vector<SatData> GNSSBatches;

    class DeltaDop
    {
        public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        unsigned int SatNum;
        Eigen::Vector3d los;
        double deltaD;
    };

    typedef std::vector<DeltaDop> DeltaDops;

    class VelocityData
    {
        public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        Eigen::Vector3d vel;
        bool valid;
    };

    class Filter
    {
        public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        Filter();
        bool Initialize( Eigen::Quaterniond& init_rot, Eigen::Matrix3d& init_rot_cov, Eigen::Vector3d& acc_pre, Eigen::Vector3d& gyr_pre );
        bool PushBack( D2R::IMUBatches& imu_batches );
        bool Update( std::vector<D2R::DeltaDop>& delta_dops, IntegrationBase* integration_ptr );

        Eigen::Quaterniond cur_rot;
        Eigen::Vector3d cur_gyr_bias;
        Eigen::Matrix<double, 6, 6> covariance;
        Eigen::Vector3d acc_pre;
        Eigen::Vector3d gyr_pre;
        double dt_pre;
        double run_time;
        bool initial_flag;
    };

    class UnscentedFilter
    {
        public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        UnscentedFilter();
        bool Initialize( Eigen::Quaterniond& init_rot, Eigen::Matrix3d& init_rot_cov, Eigen::Vector3d& acc_pre, Eigen::Vector3d& gyr_pre );
        bool PushBack( IMUBatches& imu_batches );
        bool Update( Eigen::Vector3d& delta_v_local, IntegrationBase* integration_ptr, Eigen::Matrix3d& delta_v_local_cov );
        bool Update( std::vector<DeltaDop>& delta_dops, IntegrationBase* integration_ptr );

        Eigen::Quaterniond cur_rot;
        Eigen::Vector3d cur_gyr_bias;
        Eigen::Matrix<double, 6, 6> covariance;
        Eigen::Vector3d acc_pre;
        Eigen::Vector3d gyr_pre;
        double dt_pre;
        double run_time;
        double initial_flag;
    };

    class Transformer
    {
        public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        Transformer();
        ~Transformer();
        void Initialize( IMUBatches& imu_batches, GNSSBatches& gnss_batches, unsigned int win_size, Eigen::Vector3d& usr_pos );
        void PushBack( IMUBatches& imu_batches, GNSSBatches& gnss_batches );
        bool CalcAttitude( Eigen::Matrix3d& Ceb, Eigen::Matrix3d& Cov );
        bool CalcAccuracy( Eigen::Matrix3d& Ceb, Eigen::Matrix3d& Cov );

        static void GNSSDoppler2Velocity( Eigen::Vector3d& pos, GNSSBatches& gnss_batches, Eigen::Vector4d& vel, Eigen::Matrix4d& degree_of_precision );

        GNSSBatches gnss_batches_pre;
        IMUBatches imu_batches_pre;
        IntegrationBase* integration_cur_ptr;
        bool initialize_flag;
        unsigned int win_size;
        std::queue<IntegrationBase*> integration_btw;
        std::queue<std::vector<DeltaDop>> delta_dop_btw;
        Eigen::Vector3d usr_pos;
        Eigen::Vector3d usr_pos_ecef;
        Eigen::Vector3d gravity_ecef;

        Filter state_filter;
        UnscentedFilter unscented_state_filter;

        int exec_cnt;
    };
}
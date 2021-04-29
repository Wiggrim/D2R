#include "Dop2Rot.h"

namespace D2R
{
    Transformer::Transformer()
    {
        this->initialize_flag = false;
        this->integration_cur_ptr = nullptr;
        this->win_size = 0;
        this->exec_cnt = 0;
    }

    Transformer::~Transformer()
    {
        while( !integration_btw.empty() )
        {
            IntegrationBase* integration_ptr = integration_btw.front();
            integration_btw.pop();
            delete integration_ptr;
        }
    }

    void Transformer::Initialize( IMUBatches& imu_batches, GNSSBatches& gnss_batches, unsigned int win_size, Eigen::Vector3d& usr_pos )
    {
        this->win_size = win_size;
        this->usr_pos = usr_pos;    // lat, lon, hae

        double earth_model_rn, earth_model_rm;
        double pos_ecef[3];
        double Cen[9];
        double gravity_n[3];
        double gravity_e[3];
        double local_gravity;
        Utility::CalcR( &earth_model_rm, &earth_model_rn, EM_a, EM_e2, usr_pos.x() );
        Utility::Pos2ECEF( pos_ecef, EM_e2, earth_model_rn, usr_pos.x(), usr_pos.y(), usr_pos.z() );
        this->usr_pos_ecef = Eigen::Vector3d( pos_ecef[0], pos_ecef[1], pos_ecef[2] );
        local_gravity = Utility::CalcGravity( usr_pos.x(), usr_pos.z() );
        Utility::Pos2DCM( Cen, usr_pos.x(), usr_pos.y() );
        gravity_n[0] = 0.; gravity_n[1] = 0.; gravity_n[2] = local_gravity;
        Utility::MatrixMultipy( gravity_e, Cen, 3, 3, gravity_n, 3, 1 );
        this->gravity_ecef = Eigen::Vector3d( gravity_e[0], gravity_e[1], gravity_e[2] );
        std::cout << "Derived Gravity in ECEF : \n" << this->gravity_ecef << std::endl;

        unsigned int imu_size = imu_batches.size();
        Eigen::Vector3d acc0 = imu_batches[imu_size-1].acc;
        Eigen::Vector3d gyr0 = imu_batches[imu_size-1].gyr;
        
        this->gnss_batches_pre = gnss_batches;
        this->imu_batches_pre = imu_batches;

        std::cout << "Finish Initialize with Parameter: " << std::endl;
        std::cout << "Window Size : " << win_size << std::endl;
        std::cout << "Usr Position : \n" << usr_pos << std::endl;
        std::cout << "Last IMU Measurement : \n" << acc0 << std::endl << gyr0 << std::endl;
        std::cout << "Previous GNSS Satellite Number : " << this->gnss_batches_pre.size() << std::endl;

        this->exec_cnt = 0;
        this->initialize_flag = true;

        return;
    }

    void Transformer::PushBack( IMUBatches& imu_batches, GNSSBatches& gnss_batches )
    {
        if( this->initialize_flag == false ) 
        {
            std::cout << "Need to be called after initialize\n";
            return;
        }

        unsigned int imu_size_pre = this->imu_batches_pre.size();
        Eigen::Vector3d acc0 = imu_batches_pre[imu_size_pre-1].acc;
        Eigen::Vector3d gyr0 = imu_batches_pre[imu_size_pre-1].gyr;
        this->integration_cur_ptr = new IntegrationBase( acc0, gyr0, Eigen::Vector3d(0., 0., 0.), Eigen::Vector3d( this->unscented_state_filter.cur_gyr_bias ));

        unsigned int imu_size = imu_batches.size();
        for( unsigned int i = 0; i < imu_size; i++ )
        {
            this->integration_cur_ptr->push_back( imu_batches[i].dt, imu_batches[i].acc, imu_batches[i].gyr );
        }

        if( this->integration_btw.size() == this->win_size )
        {
            delete this->integration_btw.front();
            this->integration_btw.pop();
            this->integration_btw.push( this->integration_cur_ptr );
        }
        else if( this->integration_btw.size() < this->win_size )
        {
            this->integration_btw.push( this->integration_cur_ptr );
        }
        else
        {
            std::cout << "Invalid IMU integration count : \n" << this->integration_btw.size() << ", " << this->win_size << std::endl;
            return;
        }

        std::vector<DeltaDop> dops;
        dops.clear();
        for( auto& sat_data : gnss_batches )
        {
            GNSSBatches::iterator it = std::find( gnss_batches_pre.begin(), gnss_batches_pre.end(), sat_data );
            if( it != gnss_batches_pre.end() )
            {
                Eigen::Vector3d sat_pos_pre = it->satPos;
                Eigen::Vector3d sat_pos_cur = sat_data.satPos;
                Eigen::Vector3d sat_vel_pre = it->satVel;
                Eigen::Vector3d sat_vel_cur = sat_data.satVel;
                Eigen::Vector3d los_pre = ( sat_pos_pre - usr_pos_ecef ) / (sat_pos_pre - usr_pos_ecef).norm();
                Eigen::Vector3d los_cur = ( sat_pos_cur - usr_pos_ecef ) / (sat_pos_cur - usr_pos_ecef).norm();
                Eigen::Vector3d los_average = ( los_pre + los_cur ) / 2.;
                double delta_d = it->carrierDoppler / it->carrierFrequency * PAR_C - sat_data.carrierDoppler / sat_data.carrierFrequency * PAR_C - \
                                                    los_pre.dot( it->satVel ) + los_cur.dot( sat_data.satVel ) - los_average.dot( gravity_ecef ) * this->integration_cur_ptr->sum_dt;

                DeltaDop delta_dop;
                delta_dop.deltaD = delta_d;
                delta_dop.los = los_average;
                delta_dop.SatNum = it->satNum;
                dops.push_back( delta_dop );
            }
        }

        if( delta_dop_btw.size() == this->win_size )
        {
            delta_dop_btw.pop();
            delta_dop_btw.push( dops );
        }
        else if( delta_dop_btw.size() < this->win_size )
        {
            delta_dop_btw.push( dops );
        }
        else
        {
            std::cout << "Invalid Delta Doppler Size : \n" << delta_dop_btw.size() << ", " << this->win_size << std::endl;
        }
        
        if( delta_dop_btw.size() != integration_btw.size() )
        {
            std::cout << "Miss Match between Delta Doppler and IMU Integration\n" << delta_dop_btw.size() << ", " << integration_btw.size() << std::endl;
        }

        // Change of velocity with respect to inertial frame
        Eigen::Vector4d velocity_pre, velocity_cur;
        Eigen::Matrix4d degree_of_precision;
        GNSSDoppler2Velocity( this->usr_pos, gnss_batches_pre, velocity_pre, degree_of_precision );
        GNSSDoppler2Velocity( this->usr_pos, gnss_batches, velocity_cur, degree_of_precision );
        bool velocity_update_valid;
        Eigen::Vector3d delta_v_local;
        Eigen::Matrix3d delta_v_local_cov;
        if( velocity_pre(3) != 0. && velocity_cur(3) != 0. )
        {
            velocity_update_valid = true;
            Eigen::Vector3d average_velocity = ( velocity_pre.block( 0, 0, 3, 1 ) + velocity_cur.block( 0, 0, 3, 1 ) ) / 2.;
            Eigen::Matrix3d skew_average_velocity;
            Utility::SkewMatrix( average_velocity, skew_average_velocity );
            Eigen::Vector3d earth_rot_vec = Eigen::Vector3d( 0., 0., 1. ) * EM_we * integration_btw.back()->sum_dt;
            Eigen::Quaterniond earth_rot_q = Eigen::Quaterniond( 1., earth_rot_vec(0) / 2., earth_rot_vec(1) / 2., earth_rot_vec(2) / 2. );
            earth_rot_q.normalize();
            delta_v_local = earth_rot_q.toRotationMatrix() * velocity_cur.block( 0, 0, 3, 1 ) - velocity_pre.block( 0, 0, 3, 1 ) - skew_average_velocity * earth_rot_vec - gravity_ecef * integration_btw.back()->sum_dt;
            delta_v_local_cov = degree_of_precision.block( 0, 0, 3, 3 ) * DOP_STD * DOP_STD;
        }

        if( unscented_state_filter.initial_flag == true )
        {
            // =============== UKF for mixed update ===============
            // if( velocity_update_valid == true )
            // {
            //     unscented_state_filter.Update( delta_v_local, integration_btw.back(), delta_v_local_cov );
            //     unscented_state_filter.PushBack( imu_batches );
            // }
            // else if( delta_dop_btw.back().size() != 0 )
            // {
            //     unscented_state_filter.Update( delta_dop_btw.back(), integration_btw.back() );
            //     unscented_state_filter.PushBack( imu_batches );
            // }
            // else
            // {
            //     unscented_state_filter.PushBack( imu_batches );
            // }
            // =============== UKF for Dop update ===============
            unscented_state_filter.Update( delta_dop_btw.back(), integration_btw.back() );
            unscented_state_filter.PushBack( imu_batches );            
            // =============== UKF for velocity update
            // if( velocity_update_valid == true )
            // {
            //     unscented_state_filter.Update( delta_v_local, integration_btw.back(), delta_v_local_cov );
            //     unscented_state_filter.PushBack( imu_batches );
            // }
            // else
            // {
            //     unscented_state_filter.PushBack( imu_batches );
            // }            
        }

        imu_batches_pre = imu_batches;
        gnss_batches_pre = gnss_batches;
        exec_cnt = exec_cnt + 1;

        return;
    }

    void CalcDeltaV_ECEF( std::vector<DeltaDop>& delta_dops, IntegrationBase* integration_ptr, \
                                                    std::vector<Eigen::Vector3d, Eigen::aligned_allocator<Eigen::Vector3d>>& delta_v_e, std::vector<unsigned int>& valid, unsigned int num )
    {
        unsigned int satellite_number = delta_dops.size();
        Eigen::MatrixXd los_matrix = Eigen::MatrixXd::Zero( satellite_number, 3 );
        Eigen::VectorXd delta_d_vector = Eigen::MatrixXd::Zero( satellite_number, 1 );
        for( unsigned int i = 0; i < satellite_number; i++ )
        {
            los_matrix.block( i, 0, 1, 3 ) = delta_dops[i].los.transpose();
            delta_d_vector( i, 0 ) = delta_dops[i].deltaD;
        }
        Eigen::JacobiSVD<Eigen::MatrixXd> los_matrix_svd( los_matrix, Eigen::ComputeThinU | Eigen::ComputeThinV );
        double los_matrix_cond = los_matrix_svd.singularValues()(0) / los_matrix_svd.singularValues()(los_matrix_svd.singularValues().size()-1);
        if( los_matrix_cond <= 10. && satellite_number >= 3 )
        {
            delta_v_e[num] = los_matrix_svd.solve( delta_d_vector );
            valid[num] = 1;
        }
        else
        {
            valid[num] = 0;
        }

        return;
    }

    bool Transformer::CalcAttitude(Eigen::Matrix3d& Ceb, Eigen::Matrix3d& Cov)
    {
        if( integration_btw.size() != win_size || delta_dop_btw.size() != win_size || integration_btw.size() != delta_dop_btw.size() )
        {
            std::cout << "Valid Window Length not Equal to Preset Window Size : " << integration_btw.size() << ", " << delta_dop_btw.size() << std::endl;
            return false;
        }

        bool valid_win = false;
        Eigen::Vector3d delta_v_ecef;
        Eigen::Vector3d delta_v_body;
        Eigen::Quaterniond rot_btw = Eigen::Quaterniond( 1., 0., 0., 0. );
        Eigen::Matrix3d cov = Eigen::Matrix3d::Zero();
        Eigen::Matrix3d ex_cov = Eigen::Matrix3d::Zero();

        // std::vector<Eigen::Vector3d, Eigen::aligned_allocator<Eigen::Vector3d>> delta_v_e;
        // std::vector<Eigen::Vector3d, Eigen::aligned_allocator<Eigen::Vector3d>> delta_v_b;
        // std::vector<unsigned int> valid;
        // delta_v_e.clear(); delta_v_e.resize( win_size );
        // delta_v_b.clear(); delta_v_b.resize( win_size );
        // valid.clear(); valid.resize( win_size );
        // std::vector<std::thread> threads;
        // threads.clear();

        IntegrationBase* integration_in_win_ptr = nullptr;
        std::vector<DeltaDop> delta_dop_in_win;
        double sum_dt = 0.;

        for( unsigned int i = 0; i < win_size; i++ )
        {
            valid_win = false;
            integration_in_win_ptr = integration_btw.front(); integration_btw.pop(); integration_btw.push( integration_in_win_ptr );
            delta_dop_in_win = delta_dop_btw.front(); delta_dop_btw.pop(); delta_dop_btw.push( delta_dop_in_win );
            unsigned int valid_sat_count = delta_dop_in_win.size();

            if( valid_sat_count == 0 ) continue;

            //std::thread t( &CalcDeltaV_ECEF, std::ref(delta_dop_in_win), integration_in_win_ptr, std::ref(delta_v_e), std::ref(valid), i );
            //threads.push_back( std::move(t) );

            Eigen::MatrixXd los_matrix = Eigen::MatrixXd::Zero( valid_sat_count, 3 );
            Eigen::MatrixXd delta_d_vector = Eigen::MatrixXd::Zero( valid_sat_count, 1 );
            for( unsigned int j = 0; j < valid_sat_count; j++ )
            {
                los_matrix.block( j, 0, 1, 3 ) = delta_dop_in_win[j].los.transpose();
                delta_d_vector( j, 0 ) = delta_dop_in_win[j].deltaD;
            }
            Eigen::JacobiSVD<Eigen::MatrixXd> los_matrix_svd( los_matrix, Eigen::ComputeThinU | Eigen::ComputeThinV );
            double los_matrix_cond = los_matrix_svd.singularValues()(0) / los_matrix_svd.singularValues()(los_matrix_svd.singularValues().size()-1);
            if( los_matrix_cond <= 10. && valid_sat_count >= 3 )
            {
                Eigen::Vector3d earth_rot_vec = Eigen::Vector3d( 0., 0., 1. ) * EM_we * sum_dt;
                Eigen::Quaterniond earth_rot_q = Eigen::Quaterniond( 1., earth_rot_vec(0) / 2., earth_rot_vec(1) / 2., earth_rot_vec(2) / 2. );
                earth_rot_q.normalize();
                delta_v_ecef = earth_rot_q.toRotationMatrix() * los_matrix_svd.solve( delta_d_vector );
                
                valid_win = true;
            }
            else
            {
                valid_win = false;
            }

            delta_v_body = rot_btw.toRotationMatrix() * integration_in_win_ptr->delta_v;
            rot_btw = rot_btw * integration_in_win_ptr->delta_q;

            

            if( valid_win == true )
            {
                delta_v_ecef = delta_v_ecef / delta_v_ecef.norm() * delta_v_body.norm();
                cov += delta_v_body * delta_v_body.transpose();
                ex_cov += delta_v_body * delta_v_ecef.transpose();
            }

            sum_dt += integration_in_win_ptr->sum_dt;
        }

        // for( auto& thread_t : threads )
        // {
        //     thread_t.join();
        // }

        //cov = Eigen::Matrix3d::Zero();
        //ex_cov = Eigen::Matrix3d::Zero();
        // for( unsigned int i = 0; i < win_size; i++ )
        // {
        //     if( valid[i] == 1 )
        //     {
        //         delta_v_e[i] = delta_v_e[i] / delta_v_e[i].norm() * delta_v_b[i].norm();
        //         cov += delta_v_b[i] * delta_v_b[i].transpose();
        //         ex_cov += delta_v_b[i] * delta_v_e[i].transpose();
        //     }
        // }
        
        Eigen::Vector3d earth_rot_vec = Eigen::Vector3d( 0., 0., 1. ) * EM_we * sum_dt;
        Eigen::Quaterniond earth_rot_q = Eigen::Quaterniond( 1., earth_rot_vec(0) / 2., earth_rot_vec(1) / 2., earth_rot_vec(2) / 2. );
        earth_rot_q.normalize();

        Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> cov_eigen( cov );
        double cov_eigen_check = cov_eigen.eigenvalues()(2) / cov_eigen.eigenvalues()(1);

        if( cov_eigen_check <= 912 )
        {
            Eigen::JacobiSVD<Eigen::Matrix3d> ex_cov_svd( ex_cov, Eigen::ComputeFullU | Eigen::ComputeFullV );
            Eigen::Matrix3d U = ex_cov_svd.matrixU();
            Eigen::Matrix3d V = ex_cov_svd.matrixV();
            Ceb = V*U.transpose();

            if( Ceb.determinant() < 0. )
            {
                V.block( 0, 2, 3, 1 ) = -V.block( 0, 2, 3, 1 );
                Ceb = V*U.transpose();
            }

            Ceb = earth_rot_q.toRotationMatrix().transpose() * Ceb * rot_btw.toRotationMatrix();
            CalcAccuracy( Ceb, Cov );

            if( unscented_state_filter.initial_flag == false )
            {
                Eigen::Quaterniond qeb( Ceb );
                unscented_state_filter.Initialize( qeb, Cov, imu_batches_pre.back().acc, imu_batches_pre.back().gyr );
            }
            else
            {
                Ceb = unscented_state_filter.cur_rot.toRotationMatrix();
                Cov = unscented_state_filter.covariance.block<3,3>( 0, 0 );
            }

            return true;
        }
        else
        {
            if( unscented_state_filter.initial_flag == true )
            {
                Ceb = unscented_state_filter.cur_rot.toRotationMatrix();
                Cov = unscented_state_filter.covariance.block<3,3>( 0, 0 );

                return true;
            }
            else
            {
                return false;
            }
        }
    }

    bool Transformer::CalcAccuracy( Eigen::Matrix3d& Ceb, Eigen::Matrix3d& Cov )
    {
        if( integration_btw.size() < win_size || delta_dop_btw.size() < win_size || integration_btw.size() != delta_dop_btw.size() )
        {
            std::cout << "Invalid Data Size : " << integration_btw.size() << ", " << delta_dop_btw.size() << std::endl;
            return false;
        }

        std::vector<std::vector<Eigen::Vector3d, Eigen::aligned_allocator<Eigen::Vector3d>>> los_sets;
        std::vector<std::vector<unsigned int>> sat_num_sets;
        std::vector<Eigen::Matrix3d, Eigen::aligned_allocator<Eigen::Matrix3d>> gyr_delta_v_sets;
        std::vector<Eigen::Matrix3d, Eigen::aligned_allocator<Eigen::Matrix3d>> gyr_delta_q_sets;
        std::vector<Eigen::Matrix3d, Eigen::aligned_allocator<Eigen::Matrix3d>> acc_delta_v_sets;
        std::vector<Eigen::Vector3d, Eigen::aligned_allocator<Eigen::Vector3d>> delta_v_sets;
        std::vector<Eigen::Matrix3d, Eigen::aligned_allocator<Eigen::Matrix3d>> rots;
        std::vector<Eigen::Matrix3d, Eigen::aligned_allocator<Eigen::Matrix3d>> cov_delta_q_sets;
        std::vector<Eigen::Matrix3d, Eigen::aligned_allocator<Eigen::Matrix3d>> cov_delta_v_sets;
        unsigned int total_los_num = 0;

        for( unsigned int i = 0; i < win_size; i++ )
        {
            IntegrationBase* integration_ptr = integration_btw.front(); integration_btw.pop(); integration_btw.push( integration_ptr );
            std::vector<DeltaDop> delta_dops = delta_dop_btw.front(); delta_dop_btw.pop(); delta_dop_btw.push( delta_dops );

            std::vector<Eigen::Vector3d, Eigen::aligned_allocator<Eigen::Vector3d>> los_set;
            std::vector<unsigned int> sat_num_set;

            gyr_delta_v_sets.push_back( integration_ptr->jacobian.block( 6, 12, 3, 3 ) );
            gyr_delta_q_sets.push_back( integration_ptr->jacobian.block( 3, 12, 3, 3 ) );
            acc_delta_v_sets.push_back( integration_ptr->jacobian.block( 6, 9, 3, 3 ) );
            rots.push_back( integration_ptr->delta_q.toRotationMatrix() );
            delta_v_sets.push_back( integration_ptr->delta_v );
            cov_delta_q_sets.push_back( integration_ptr->covariance.block( 3, 3, 3, 3 ) );
            cov_delta_v_sets.push_back( integration_ptr->covariance.block( 6, 6, 3, 3 ) );

            for( DeltaDop& delta_dop : delta_dops )
            {
                los_set.push_back( delta_dop.los );
                sat_num_set.push_back( delta_dop.SatNum );
                total_los_num++;
            }

            los_sets.push_back( los_set );
            sat_num_sets.push_back( sat_num_set );
        }

        Eigen::MatrixXd gyr_concate_delta_v = Eigen::MatrixXd::Zero( 3*win_size, 3 );
        Eigen::MatrixXd acc_concate_delta_v = Eigen::MatrixXd::Zero( 3*win_size, 3 );
        Eigen::MatrixXd cov_concate_delta_v = Eigen::MatrixXd::Zero( 3*win_size, 3*win_size );
        Eigen::MatrixXd converse = Eigen::MatrixXd::Zero( 3*win_size, 3*win_size );
        Eigen::Matrix3d total_rot = Eigen::Matrix3d::Identity();
        Eigen::Matrix3d gyr_total_rot = Eigen::Matrix3d::Zero();
        Eigen::Matrix3d cov_total_rot = Eigen::Matrix3d::Zero();
        std::vector<Eigen::Matrix3d, Eigen::aligned_allocator<Eigen::Matrix3d>> cov_total_rot_sets;
        std::vector<Eigen::Matrix3d, Eigen::aligned_allocator<Eigen::Matrix3d>> skew_local_delta_v_sets;
        cov_total_rot_sets.resize( win_size );
        skew_local_delta_v_sets.resize( win_size );
        Eigen::MatrixXd H = Eigen::MatrixXd::Zero( 3*win_size, 3 );
        for( int i = win_size-1; i >= 0; i-- )
        {
            gyr_total_rot += total_rot.transpose() * gyr_delta_q_sets[i];
            cov_total_rot += total_rot.transpose() * cov_delta_q_sets[i] * total_rot;
            cov_total_rot_sets[i] = cov_total_rot;
            total_rot = rots[i] * total_rot;

            Eigen::Vector3d local_delta_v = total_rot.transpose() * delta_v_sets[i];
            Eigen::Matrix3d skew_local_delta_v;
            skew_local_delta_v << 0., -local_delta_v(2), local_delta_v(1), \
                                                            local_delta_v(2), 0., -local_delta_v(0), \
                                                            -local_delta_v(1), local_delta_v(0), 0.;
            skew_local_delta_v_sets[i] = skew_local_delta_v;
            
            gyr_concate_delta_v.block( 3*i, 0, 3, 3 ) = skew_local_delta_v * gyr_total_rot + total_rot.transpose() * gyr_delta_v_sets[i];
            acc_concate_delta_v.block( 3*i, 0, 3, 3 ) = total_rot.transpose() * acc_delta_v_sets[i];
            cov_concate_delta_v.block( 3*i, 3*i, 3, 3 ) = skew_local_delta_v * cov_total_rot * skew_local_delta_v.transpose() + total_rot.transpose() * cov_delta_v_sets[i] * total_rot;
            for( unsigned int j = i+1; j < win_size; j++ )
            {
                cov_concate_delta_v.block( 3*i, 3*j, 3, 3 ) = skew_local_delta_v * cov_total_rot_sets[j] * skew_local_delta_v_sets[j].transpose();
                cov_concate_delta_v.block( 3*j, 3*i, 3, 3 ) = skew_local_delta_v_sets[j].transpose() * cov_total_rot_sets[j] * skew_local_delta_v.transpose();
            }
            converse.block( 3*i, 3*i, 3, 3 ) = Ceb;
            H.block( 3*i, 0, 3, 3 ) = Ceb * skew_local_delta_v;
        }

        Eigen::MatrixXd cov_delta_d = Eigen::MatrixXd::Zero( total_los_num, total_los_num );
        Eigen::MatrixXd delta_v_delta_d = Eigen::MatrixXd::Zero( 3*win_size, total_los_num );
        std::vector<unsigned int> valid_win_index;
        unsigned int valid_delta_v_size = 0;
        unsigned int passed_num = 0;
        for( unsigned int i = 0; i < win_size; i++ )
        {
            std::vector<Eigen::Vector3d, Eigen::aligned_allocator<Eigen::Vector3d>> los_set;
            std::vector<unsigned int> sat_num_set;

            los_set = los_sets[i];
            sat_num_set = sat_num_sets[i];

            unsigned int valid_num = los_set.size();
            if( valid_num == 0 ) continue;

            Eigen::MatrixXd los_mat = Eigen::MatrixXd::Zero( valid_num, 3 );
            for( unsigned int j = 0; j < valid_num; j++ )
            {
                los_mat.block( j, 0, 1, 3 ) = los_set[j].transpose();
            }
            cov_delta_d.block( passed_num, passed_num, valid_num, valid_num ) = Eigen::MatrixXd::Identity( valid_num, valid_num ) * 2 * DOP_STD * DOP_STD + Eigen::MatrixXd::Ones( valid_num, valid_num ) * OSC_STD * OSC_STD;

            if( i != win_size-1 )
            {
                std::vector<unsigned int> sat_num_set_next;
                sat_num_set_next = sat_num_sets[i+1];
                for( unsigned int j = 0; j < valid_num; j++ )
                {
                    unsigned int sat_num = sat_num_set[j];
                    std::vector<unsigned int>::iterator it = std::find( sat_num_set_next.begin(), sat_num_set_next.end(), sat_num );
                    if( it != sat_num_set_next.end() )
                    {
                        unsigned int next_index = it - sat_num_set_next.begin();
                        cov_delta_d( passed_num+j, passed_num+valid_num+next_index ) = - DOP_STD * DOP_STD;
                        cov_delta_d( passed_num+valid_num+next_index, passed_num+j ) = - DOP_STD * DOP_STD;
                    }
                }
            }

            Eigen::JacobiSVD<Eigen::MatrixXd> los_svd( los_mat );
            double los_cond_num = los_svd.singularValues()(0) / los_svd.singularValues()(los_svd.singularValues().size()-1);
            if( los_cond_num <= 10. && los_mat.rows() >= 3 )
            {
                delta_v_delta_d.block( valid_delta_v_size, passed_num, 3, valid_num ) = ( los_mat.transpose() * los_mat ).inverse() * los_mat.transpose();
                valid_delta_v_size += 3;
                valid_win_index.push_back( i );
            }
            passed_num += valid_num;
        }

        unsigned int valid_win_num = valid_win_index.size();
        Eigen::MatrixXd sort_mat = Eigen::MatrixXd::Zero( valid_win_num*3, win_size*3 );
        for( unsigned int i = 0; i < valid_win_num; i++ )
        {
            unsigned int win_index = valid_win_index[i];
            sort_mat.block( 3*i, 3*win_index, 3, 3 ) = Eigen::Matrix3d::Identity();
        }


        Eigen::MatrixXd cov_in = delta_v_delta_d.block(0, 0, valid_win_num*3, total_los_num) * cov_delta_d * delta_v_delta_d.block(0, 0, valid_win_num*3, total_los_num).transpose() + \
                                                    sort_mat * converse * ( cov_concate_delta_v + gyr_concate_delta_v*Eigen::Matrix3d::Identity()*gyr_concate_delta_v.transpose()*INIT_GYR_STD * INIT_GYR_STD + \
                                                    acc_concate_delta_v*Eigen::Matrix3d::Identity()*acc_concate_delta_v.transpose()*INIT_ACC_STD * INIT_ACC_STD ) * converse.transpose() * sort_mat.transpose();
        Eigen::MatrixXd H_inverse = ( ( sort_mat * H ).transpose() * sort_mat * H ).inverse() * (sort_mat * H).transpose();
        
        Cov = H_inverse * cov_in * H_inverse.transpose();

        return true;
    }

    void Transformer::GNSSDoppler2Velocity( Eigen::Vector3d& pos, GNSSBatches& gnss_batches, Eigen::Vector4d& vel, Eigen::Matrix4d& degree_of_precision )
    {
        unsigned int sat_meas_num = gnss_batches.size();
        if( sat_meas_num < 4 )
        {
            std::cout << "Less than 4 satellites for velocity derivation\n";
            vel.setZero();
            return;
        }

        Eigen::Vector3d pos_ecef;
        double Rm, Rn;
        Utility::CalcR( &Rm, &Rn, EM_a, EM_e2, pos.x() );
        Utility::Pos2ECEF( pos_ecef.data(), EM_e2, Rn, pos.x(), pos.y(), pos.z() );

        Eigen::MatrixXd y = Eigen::MatrixXd::Zero( sat_meas_num, 1 );
        Eigen::MatrixXd H = Eigen::MatrixXd::Zero( sat_meas_num, 4 );
        for( unsigned int i = 0; i < sat_meas_num; i++ )
        {
            SatData& sat_data = gnss_batches[i];
            Eigen::Vector3d los = ( sat_data.satPos - pos_ecef ).normalized();
            H.block( i, 0, 1, 3 ) = -los.transpose(); H( i, 3 ) = 1.;
            y( i, 0 ) = sat_data.carrierDoppler / sat_data.carrierFrequency * PAR_C - los.dot( sat_data.satVel );
        }

        degree_of_precision = ( H.transpose() * H ).inverse();

        Eigen::JacobiSVD<Eigen::MatrixXd> H_svd( H, Eigen::ComputeFullV | Eigen::ComputeFullU );
        vel = H_svd.solve( y );

        return;
    }
}
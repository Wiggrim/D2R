#include <eigen3/Eigen/Dense>
#include <vector>
#include <queue>
#include <algorithm>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <chrono>
#include <stdlib.h>
#include "Dop2Rot.h"

unsigned int used_satellite_num[3] = { 2, 5, 12 };

// =============== Start Preset vehicle's position to surrounding cities ===============
// 114.939 	40.764 Zhang jia kou
// 118.026 	41.033 Cheng de
// 116.847 	38.274 Cang zhou
// 0.698156713 2.030342625 105.906359 Actual position
Eigen::Vector3d preset_pos( 0.698156713, 2.030342625, 105.906359 );
// =============== End Preset vehicle's position to surrounding cities ===============

bool GetDataFromLine( std::string& line, std::vector<double>& data_out, unsigned int data_num, char seperate )
{
    data_out.clear();
    data_out.resize( data_num );
    std::stringstream ss( line );
    std::string subline;
    for( unsigned int i = 0; i < data_num; i++ )
    {
        if( ss.eof() )
        {
            std::cout << "Unexpected End in a line : " << line << std::endl;
            return false;
        }        
        std::getline( ss, subline, seperate );
        data_out[i] = atof( subline.c_str() );
    }

    return true;
}

int main( int argc, char** argv )
{
    if( argc != 4 ) 
    {
        std::cout << "1. Path to DataSet, 2. Window Size, 3. Output Path\n";
        return 0;
    }

    double Cen[9];
    Utility::Pos2DCM( Cen, preset_pos.x(), preset_pos.y() );
    Eigen::Matrix3d Cen_m;
    Cen_m << Cen[0], Cen[1], Cen[2], \
                        Cen[3], Cen[4], Cen[5], \
                        Cen[6], Cen[7], Cen[8];

    D2R::Transformer transformer;
    Eigen::Matrix3d Ceb, Cov;

    std::string data_path = argv[1];
    std::string output_path = argv[3];
    unsigned int win_size = atoi( argv[2] );

    D2R::IMUBatches imu_batches;
    D2R::GNSSBatches gnss_batches;
    imu_batches.clear();
    gnss_batches.clear();

    std::fstream fout( output_path, std::ios::out );
    std::fstream fin( data_path, std::ios::in );
    std::string line;
    std::string subline;
    std::vector<double> data_loads;
    data_loads.clear();
    while( !fin.eof() )
    {
        std::getline( fin, line );
        if( line.substr( 0, 4 ) == "$imu" )
        {
            subline = line.substr( 5 );
            GetDataFromLine( subline, data_loads, 6, ' ' );
            D2R::IMUMeas imu_meas;
            imu_meas.acc = Eigen::Vector3d( data_loads[0], data_loads[1], data_loads[2] );
            imu_meas.gyr = Eigen::Vector3d( data_loads[3], data_loads[4], data_loads[5] );
            imu_meas.dt = 0.01;
            imu_batches.push_back( imu_meas );
        }
        else if( line.substr( 0, 4 ) == "$gps" )
        {
            subline = line.substr( 5 );
            GetDataFromLine( subline, data_loads, 9, ' ' );
            Eigen::Vector3d pos( data_loads[0], data_loads[1], data_loads[2] );
            Eigen::Vector3d vel( data_loads[3], data_loads[4], data_loads[5] );
            double tow = data_loads[6];
            unsigned int satellite_number = data_loads[8];
            gnss_batches.resize( satellite_number );
            for( unsigned int i = 0; i < satellite_number; i++ )
            {
                std::getline( fin, line );
                GetDataFromLine( line, data_loads, 3, ' ' );
                gnss_batches[i].satPos = Eigen::Vector3d( data_loads[0], data_loads[1], data_loads[2] );
            }
            for( unsigned int i = 0; i < satellite_number; i++ )
            {
                std::getline( fin, line );
                GetDataFromLine( line, data_loads, 3, ' ' );
                gnss_batches[i].satVel = Eigen::Vector3d( data_loads[0], data_loads[1], data_loads[2] );
            }
            std::getline( fin, line );
            GetDataFromLine( line, data_loads, satellite_number, ' ' );
            for( unsigned int i = 0; i < satellite_number; i++ )
            {
                gnss_batches[i].transmitTime = data_loads[i];
            }
            std::getline( fin, line );
            GetDataFromLine( line, data_loads, satellite_number, ' ' );
            for( unsigned int i = 0; i < satellite_number; i++ )
            {
                gnss_batches[i].carrierDoppler = INVERSE_DOP ? -data_loads[i] : data_loads[i];
                gnss_batches[i].carrierFrequency = CARRIER_FREQ;
            }
            std::getline( fin, line );
            GetDataFromLine( line, data_loads, satellite_number, ' ' );
            for( unsigned int i = 0; i < satellite_number; i++ )
            {
                gnss_batches[i].satNum = data_loads[i];
            }

            std::cout << "=============== Start Attitude Estimation ===============\n";

            Eigen::Vector4d velocity_from_doppler;
            Eigen::Matrix4d degree_of_precision;
            D2R::Transformer::GNSSDoppler2Velocity( pos, gnss_batches, velocity_from_doppler, degree_of_precision );
            velocity_from_doppler.block( 0, 0, 3, 1 ) = Cen_m.transpose() * velocity_from_doppler.block( 0, 0, 3, 1 );
            std::cout << "Velocity from GNSS : " << velocity_from_doppler.transpose() << std::endl;

            //=============== Start Sorting the chosen satellites ===============
            // if( /*transformer.unscented_state_filter.run_time >= 15. && transformer.unscented_state_filter.run_time <= 75.*/ 1 )
            // {
            //     D2R::GNSSBatches::iterator it = gnss_batches.begin();
            //     while( it != gnss_batches.end() )
            //     {
            //         if( transformer.state_filter.initial_flag == false )
            //         {
            //             if( it->satNum != used_satellite_num[0] && it->satNum != used_satellite_num[1] && it->satNum != used_satellite_num[2] )
            //             {
            //                 it  = gnss_batches.erase(it);
            //             }
            //             else
            //             {
            //                 it++;
            //             }        
            //         }
            //         else
            //         {
            //             if( it->satNum != used_satellite_num[0] && it->satNum != used_satellite_num[1] && it->satNum != used_satellite_num[2]  )
            //             {
            //                 it  = gnss_batches.erase(it);
            //             }
            //             else
            //             {
            //                 it++;
            //             }    
            //         }
            //     }
            // }
            //=============== End Sorting the chosen satellites ===============

            if( transformer.initialize_flag == false )
            {
                transformer.Initialize( imu_batches, gnss_batches, win_size, preset_pos );
            }
            else
            {
                std::chrono::time_point<std::chrono::system_clock> start = std::chrono::system_clock::now();
                transformer.PushBack( imu_batches, gnss_batches );

                if( transformer.CalcAttitude( Ceb, Cov ) )
                {
                    Eigen::Matrix3d Cnb = Cen_m.transpose() * Ceb;
                    Eigen::Vector3d euler_angle;
                    Utility::Dcm2Euler( Cnb, euler_angle );
                    std::cout <<  "TOW : " << tow << std::endl << "Euler angle : " << euler_angle.transpose() \
                                                << std::endl << "Covariance : " << Cov.diagonal().transpose() << std::endl;         
                    std::cout << "Gyr Bias : " << transformer.unscented_state_filter.cur_gyr_bias.transpose() << std::endl;          
                    fout.precision(15);
                    fout << tow << ", " << euler_angle.x() << ", " << euler_angle.y() << ", " << euler_angle.z() \
                                                << ", " << Cov(0,0) << ", " << Cov(1,1) << ", " << Cov(2,2) << std::endl;
                }
                std::chrono::time_point<std::chrono::system_clock> end = std::chrono::system_clock::now();
                std::chrono::duration<double> used_seconds = end - start;
                std::cout << "Totally use " << used_seconds.count() * 1e3 << " ms" << std::endl;
            }

            imu_batches.clear();
            gnss_batches.clear();
            
        }

    }
}
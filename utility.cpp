#include "utility.h"

namespace Utility
{
	double CalcGravity(double latitude, double hae)
	{
		double a1,a2,a3,a4,a5,a6,s2,s4;
		a1 = 9.7803267714;
		a2 = 0.0052790414;
		a3 = 0.0000232718;
		a4 = -0.0000030876910891;
		a5 = 0.0000000043977311;
		a6 = 0.0000000000007211;
		s2 = sin(latitude)*sin(latitude);
		s4 = s2 * s2;
		return a1 * (1 + a2*s2 + a3*s4) + (a4 + a5*s2)*hae + a6 * hae * hae;
	}

	bool MatrixMultipy(double *c, double *a, int ar, int ac, double *b, int br, int bc)
	{
		int i,j,k;
		double t = 0.0;
		if (ac!=br)
		{
			return false;
		}
		for(i=0; i<ar; i++)
		{
			for(j=0; j<bc; j++)
			{
				for(k=0; k<ac; k++)
				{
					t = t + a[i*ac+k]*b[k*bc+j];
				}
				c[i*bc+j] = t;
				t = 0.0;
			}
		}
		return true;
	}

	void CalcR(double *Rm, double *Rn, double a, double e2, double lat)
	{
		Rm[0] = a*(1-e2)/pow((1-e2 * sin(lat)*sin(lat)),(3.0/2.0));
		Rn[0] = a/sqrt(1-e2 * sin(lat)*sin(lat));
	}

	void Pos2ECEF(double *ecef_xyz, double e2, double Rn, double lat, double lon, double hae)
	{
		double c_lat, s_lat, c_lon, s_lon, Rn_h; 
		c_lat = cos(lat);
		s_lat = sin(lat);
		c_lon = cos(lon);
		s_lon = sin(lon);

		Rn_h = Rn + hae;

		ecef_xyz[0] = Rn_h*c_lat*c_lon;
		ecef_xyz[1] = Rn_h * c_lat * s_lon;
		ecef_xyz[2] = (Rn*(1-e2)+hae)*s_lat;
	}

	void Pos2DCM(double *C_ne, double lat, double lon)
	{
		double c_lat, s_lat, c_lon, s_lon; 
		c_lat = cos(lat);
		s_lat = sin(lat);
		c_lon = cos(lon);
		s_lon = sin(lon);

		C_ne[0] = -s_lat*c_lon;
		C_ne[1] = -s_lon;
		C_ne[2] = -c_lat*c_lon;
		C_ne[3] = -s_lat*s_lon;
		C_ne[4] = c_lon;
		C_ne[5] = -c_lat*s_lon;
		C_ne[6] = c_lat;
		C_ne[7] = 0.0;
		C_ne[8] = -s_lat;
	}

	bool Dcm2Euler( Eigen::Matrix3d& C, Eigen::Vector3d& Euler )
	{
		double roll, pitch, yaw;
		pitch = atan( -C(2,0) / sqrt( C(2,1)*C(2,1) + C(2,2)*C(2,2) ) );
		roll = atan2( C(2,1), C(2,2) );
		yaw = atan2( C(1,0), C(0,0) );

		Euler = Eigen::Vector3d( roll, pitch, yaw );
		return true;
	}

	void SkewMatrix( Eigen::Vector3d& vec, Eigen::Matrix3d& skew_mat )
	{
		skew_mat << 0., -vec(2), vec(1), \
									vec(2), 0., -vec(0), \
									-vec(1), vec(0), 0.;

		return;
	}

}
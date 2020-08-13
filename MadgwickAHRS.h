//=====================================================================================================
// MadgwickAHRS.h
//=====================================================================================================
//
// Implementation of Madgwick's IMU and AHRS algorithms.
// See: http://www.x-io.co.uk/node/8#open_source_ahrs_and_imu_algorithms
//
// Date			Author			Notes
// 29/09/2011	SOH Madgwick	Initial release
// 02/10/2011	SOH Madgwick	Optimised for reduced CPU load
// 07/07/2019	Lukas			Adaptations for own projects
//=====================================================================================================
#ifndef MadgwickAHRS_h
#define MadgwickAHRS_h

class MADGWICK_AHRS {
	
	public:
	// MADGWICK_AHRS constructor
	MADGWICK_AHRS(float beta, float zeta = 0);
	
	// MADGWICK_AHRS destructor
	~MADGWICK_AHRS(void);
	
	// set beta value
	void set_beta(float beta);
	
	// set zeta value
	void set_zeta(float zeta);
	
	// get pose in euler angles
	void get_euler(float dt_s, float ax, float ay, float az, float gx, float gy, float gz, float mx, float my, float mz, float &angle_x, float &angle_y, float &angle_z);
	
	// get gyro bias estimation
	void get_gyroBias(float &gbx, float &gby, float &gbz);
	
	private:
	float m_beta;	// algorithm gain (2 * proportional gain (Kp))	(0.041 MARG, 0.033 IMU)
	float m_zeta;	// gain to remove gyroscope bias errors 		(0.003 - 0.015)
	
	float m_q0,	m_q1, m_q2, m_q3;		// quaternion of sensor frame relative to auxiliary frame
	
	// imu update time in seconds
	float m_dt_s;
	// imu measurements
	float m_ax, m_ay, m_az, m_gx, m_gy, m_gz, m_mx, m_my, m_mz;
	
	// estimated gyro biases
	float m_gbx = 0;
	float m_gby = 0;
	float m_gbz = 0;
	
	// AHRS algorithm update
	void madgwickAHRSupdate();
	// IMU algorithm update
	void madgwickAHRSupdateIMU();
	// Fast inverse square-root
	float invSqrt(float x);
};

#endif
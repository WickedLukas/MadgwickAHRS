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
	MADGWICK_AHRS(float beta);
	
	// MADGWICK_AHRS destructor
	~MADGWICK_AHRS(void);
	
	// set beta value
	void set_beta(float beta);
	
	// set z-axis angle and enable z-axis rotation
	void set_angle_z(float angle_z_new);
	
	// enable or disable z-axis rotation
	void enable_z_rotation(bool z_rotation);
	
	// get pose in euler angles
	void get_euler(float &angle_x, float &angle_y, float &angle_z, float dt_s, float ax, float ay, float az, float gx, float gy, float gz, float mx, float my, float mz);
	
private:
	
	float m_beta;	// algorithm gain (2 * proportional gain (Kp))
	
	bool m_z_rotation;	// flag for z-axis rotation
	
	float m_q0,	m_q1, m_q2, m_q3;		// quaternion of sensor frame relative to auxiliary frame
	float m_qz0, m_qz3;					// rotation quaternion for z/yaw-angle
	float m_qr0, m_qr1, m_qr2, m_qr3;	// m_q rotated with m_qz
	
	// imu update time in seconds
	float m_dt_s;
	// imu measurements
	float m_ax, m_ay, m_az, m_gx, m_gy, m_gz, m_mx, m_my, m_mz;
	
	// AHRS algorithm update
	void madgwickAHRSupdate();
	// IMU algorithm update
	void madgwickAHRSupdateIMU();
	// Fast inverse square-root
	float invSqrt(float x);
};

#endif
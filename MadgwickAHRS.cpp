//=====================================================================================================
// MadgwickAHRS.cpp
//=====================================================================================================
//
// Implementation of Madgwick's IMU and AHRS algorithms.
// See: http://www.x-io.co.uk/node/8#open_source_ahrs_and_imu_algorithms
//
// Date			Author			Notes
// 29/09/2011	SOH Madgwick	Initial release
// 02/10/2011	SOH Madgwick	Optimised for reduced CPU load
// 19/02/2012	SOH Madgwick	Magnetometer measurement is normalized
// 07/07/2019	Lukas			Adaptations for own projects
//=====================================================================================================

//---------------------------------------------------------------------------------------------------
// Header files

#include "MadgwickAHRS.h"
#include <Arduino.h>

// factor for converting a radian number to an equivalent number in degrees
const float RAD2DEG = 4068 / 71;

// MADGWICK_AHRS constructor
MADGWICK_AHRS::MADGWICK_AHRS(float beta) {
	m_beta = beta;
	m_z_rotation = false;
	
	m_q0 = 1; m_q1 = 0; m_q2 = 0; m_q3 = 0;
	m_qz0 = 1; m_qz3 = 0;
}

// MADGWICK_AHRS destructor
MADGWICK_AHRS::~MADGWICK_AHRS(void) {
}

// set beta value
void MADGWICK_AHRS::set_beta(float beta) {
	m_beta = beta;
}

// get pose in euler angles
void MADGWICK_AHRS::get_euler(float &angle_x, float &angle_y, float &angle_z, float dt_s, float ax, float ay, float az, float gx, float gy, float gz, float mx, float my, float mz) {
	m_dt_s = dt_s;
	
	m_ax = ax; m_ay = ay; m_az = az;
	m_gx = gx; m_gy = gy; m_gz = gz;
	m_mx = mx; m_my = my; m_mz = mz;
	
	madgwickAHRSupdate();
	
	// euler angles for left handed coordinate system with zyx rotation sequence
	angle_z = atan2(2*(m_q1*m_q2 + m_q0*m_q3), m_q0*m_q0 + m_q1*m_q1 - m_q2*m_q2 - m_q3*m_q3) * RAD2DEG;
	angle_y = asin(-2*(m_q1*m_q3 - m_q0*m_q2)) * RAD2DEG;
	angle_x = atan2(2*(m_q2*m_q3 + m_q0*m_q1), m_q0*m_q0 - m_q1*m_q1 - m_q2*m_q2 + m_q3*m_q3) * RAD2DEG;
}

// AHRS algorithm update
void MADGWICK_AHRS::madgwickAHRSupdate() {
	static float recipNorm;
	static float s0, s1, s2, s3;
	static float qDot1, qDot2, qDot3, qDot4;
	static float hx, hy;
	static float _2q0mx, _2q0my, _2q0mz, _2q1mx, _2bx, _2bz, _4bx, _4bz, _2q0, _2q1, _2q2, _2q3, _2q0q2, _2q2q3, q0q0, q0q1, q0q2, q0q3, q1q1, q1q2, q1q3, q2q2, q2q3, q3q3;

	// Use IMU algorithm if magnetometer measurement invalid (avoids NaN in magnetometer normalization)
	if((m_mx == 0.0f) && (m_my == 0.0f) && (m_mz == 0.0f)) {
		madgwickAHRSupdateIMU();
		return;
	}

	// Rate of change of quaternion from gyroscope
	qDot1 = 0.5f * (-m_q1 * m_gx - m_q2 * m_gy - m_q3 * m_gz);
	qDot2 = 0.5f * (m_q0 * m_gx + m_q2 * m_gz - m_q3 * m_gy);
	qDot3 = 0.5f * (m_q0 * m_gy - m_q1 * m_gz + m_q3 * m_gx);
	qDot4 = 0.5f * (m_q0 * m_gz + m_q1 * m_gy - m_q2 * m_gx);

	// Compute feedback only if accelerometer measurement valid (avoids NaN in accelerometer normalization)
	if(!((m_ax == 0.0f) && (m_ay == 0.0f) && (m_az == 0.0f))) {

		// normalize accelerometer measurement
		recipNorm = invSqrt(m_ax * m_ax + m_ay * m_ay + m_az * m_az);
		m_ax *= recipNorm;
		m_ay *= recipNorm;
		m_az *= recipNorm;

		// normalize magnetometer measurement
		recipNorm = invSqrt(m_mx * m_mx + m_my * m_my + m_mz * m_mz);
		m_mx *= recipNorm;
		m_my *= recipNorm;
		m_mz *= recipNorm;

		// Auxiliary variables to avoid repeated arithmetic
		_2q0mx = 2.0f * m_q0 * m_mx;
		_2q0my = 2.0f * m_q0 * m_my;
		_2q0mz = 2.0f * m_q0 * m_mz;
		_2q1mx = 2.0f * m_q1 * m_mx;
		_2q0 = 2.0f * m_q0;
		_2q1 = 2.0f * m_q1;
		_2q2 = 2.0f * m_q2;
		_2q3 = 2.0f * m_q3;
		_2q0q2 = 2.0f * m_q0 * m_q2;
		_2q2q3 = 2.0f * m_q2 * m_q3;
		q0q0 = m_q0 * m_q0;
		q0q1 = m_q0 * m_q1;
		q0q2 = m_q0 * m_q2;
		q0q3 = m_q0 * m_q3;
		q1q1 = m_q1 * m_q1;
		q1q2 = m_q1 * m_q2;
		q1q3 = m_q1 * m_q3;
		q2q2 = m_q2 * m_q2;
		q2q3 = m_q2 * m_q3;
		q3q3 = m_q3 * m_q3;

		// Reference direction of Earth's magnetic field
		hx = m_mx * q0q0 - _2q0my * m_q3 + _2q0mz * m_q2 + m_mx * q1q1 + _2q1 * m_my * m_q2 + _2q1 * m_mz * m_q3 - m_mx * q2q2 - m_mx * q3q3;
		hy = _2q0mx * m_q3 + m_my * q0q0 - _2q0mz * m_q1 + _2q1mx * m_q2 - m_my * q1q1 + m_my * q2q2 + _2q2 * m_mz * m_q3 - m_my * q3q3;
		_2bx = sqrt(hx * hx + hy * hy);
		_2bz = -_2q0mx * m_q2 + _2q0my * m_q1 + m_mz * q0q0 + _2q1mx * m_q3 - m_mz * q1q1 + _2q2 * m_my * m_q3 - m_mz * q2q2 + m_mz * q3q3;
		_4bx = 2.0f * _2bx;
		_4bz = 2.0f * _2bz;

		// Gradient decent algorithm corrective step
		s0 = -_2q2 * (2.0f * q1q3 - _2q0q2 - m_ax) + _2q1 * (2.0f * q0q1 + _2q2q3 - m_ay) - _2bz * m_q2 * (_2bx * (0.5f - q2q2 - q3q3) + _2bz * (q1q3 - q0q2) - m_mx) + (-_2bx * m_q3 + _2bz * m_q1) * (_2bx * (q1q2 - q0q3) + _2bz * (q0q1 + q2q3) - m_my) + _2bx * m_q2 * (_2bx * (q0q2 + q1q3) + _2bz * (0.5f - q1q1 - q2q2) - m_mz);
		s1 = _2q3 * (2.0f * q1q3 - _2q0q2 - m_ax) + _2q0 * (2.0f * q0q1 + _2q2q3 - m_ay) - 4.0f * m_q1 * (1 - 2.0f * q1q1 - 2.0f * q2q2 - m_az) + _2bz * m_q3 * (_2bx * (0.5f - q2q2 - q3q3) + _2bz * (q1q3 - q0q2) - m_mx) + (_2bx * m_q2 + _2bz * m_q0) * (_2bx * (q1q2 - q0q3) + _2bz * (q0q1 + q2q3) - m_my) + (_2bx * m_q3 - _4bz * m_q1) * (_2bx * (q0q2 + q1q3) + _2bz * (0.5f - q1q1 - q2q2) - m_mz);
		s2 = -_2q0 * (2.0f * q1q3 - _2q0q2 - m_ax) + _2q3 * (2.0f * q0q1 + _2q2q3 - m_ay) - 4.0f * m_q2 * (1 - 2.0f * q1q1 - 2.0f * q2q2 - m_az) + (-_4bx * m_q2 - _2bz * m_q0) * (_2bx * (0.5f - q2q2 - q3q3) + _2bz * (q1q3 - q0q2) - m_mx) + (_2bx * m_q1 + _2bz * m_q3) * (_2bx * (q1q2 - q0q3) + _2bz * (q0q1 + q2q3) - m_my) + (_2bx * m_q0 - _4bz * m_q2) * (_2bx * (q0q2 + q1q3) + _2bz * (0.5f - q1q1 - q2q2) - m_mz);
		s3 = _2q1 * (2.0f * q1q3 - _2q0q2 - m_ax) + _2q2 * (2.0f * q0q1 + _2q2q3 - m_ay) + (-_4bx * m_q3 + _2bz * m_q1) * (_2bx * (0.5f - q2q2 - q3q3) + _2bz * (q1q3 - q0q2) - m_mx) + (-_2bx * m_q0 + _2bz * m_q2) * (_2bx * (q1q2 - q0q3) + _2bz * (q0q1 + q2q3) - m_my) + _2bx * m_q1 * (_2bx * (q0q2 + q1q3) + _2bz * (0.5f - q1q1 - q2q2) - m_mz);
		recipNorm = invSqrt(s0 * s0 + s1 * s1 + s2 * s2 + s3 * s3); // normalize step magnitude
		s0 *= recipNorm;
		s1 *= recipNorm;
		s2 *= recipNorm;
		s3 *= recipNorm;

		// Apply feedback step
		qDot1 -= m_beta * s0;
		qDot2 -= m_beta * s1;
		qDot3 -= m_beta * s2;
		qDot4 -= m_beta * s3;
	}

	// Integrate rate of change of quaternion to yield quaternion
	m_q0 += qDot1 * m_dt_s;
	m_q1 += qDot2 * m_dt_s;
	m_q2 += qDot3 * m_dt_s;
	m_q3 += qDot4 * m_dt_s;

	// normalize quaternion
	recipNorm = invSqrt(m_q0 * m_q0 + m_q1 * m_q1 + m_q2 * m_q2 + m_q3 * m_q3);
	m_q0 *= recipNorm;
	m_q1 *= recipNorm;
	m_q2 *= recipNorm;
	m_q3 *= recipNorm;
}

// IMU algorithm update
void MADGWICK_AHRS::madgwickAHRSupdateIMU() {
	static float recipNorm;
	static float s0, s1, s2, s3;
	static float qDot1, qDot2, qDot3, qDot4;
	static float _2q0, _2q1, _2q2, _2q3, _4q0, _4q1, _4q2 ,_8q1, _8q2, q0q0, q1q1, q2q2, q3q3;

	// Rate of change of quaternion from gyroscope
	qDot1 = 0.5f * (-m_q1 * m_gx - m_q2 * m_gy - m_q3 * m_gz);
	qDot2 = 0.5f * (m_q0 * m_gx + m_q2 * m_gz - m_q3 * m_gy);
	qDot3 = 0.5f * (m_q0 * m_gy - m_q1 * m_gz + m_q3 * m_gx);
	qDot4 = 0.5f * (m_q0 * m_gz + m_q1 * m_gy - m_q2 * m_gx);

	// Compute feedback only if accelerometer measurement valid (avoids NaN in accelerometer normalization)
	if(!((m_ax == 0.0f) && (m_ay == 0.0f) && (m_az == 0.0f))) {

		// normalize accelerometer measurement
		recipNorm = invSqrt(m_ax * m_ax + m_ay * m_ay + m_az * m_az);
		m_ax *= recipNorm;
		m_ay *= recipNorm;
		m_az *= recipNorm;

		// Auxiliary variables to avoid repeated arithmetic
		_2q0 = 2.0f * m_q0;
		_2q1 = 2.0f * m_q1;
		_2q2 = 2.0f * m_q2;
		_2q3 = 2.0f * m_q3;
		_4q0 = 4.0f * m_q0;
		_4q1 = 4.0f * m_q1;
		_4q2 = 4.0f * m_q2;
		_8q1 = 8.0f * m_q1;
		_8q2 = 8.0f * m_q2;
		q0q0 = m_q0 * m_q0;
		q1q1 = m_q1 * m_q1;
		q2q2 = m_q2 * m_q2;
		q3q3 = m_q3 * m_q3;

		// Gradient decent algorithm corrective step
		s0 = _4q0 * q2q2 + _2q2 * m_ax + _4q0 * q1q1 - _2q1 * m_ay;
		s1 = _4q1 * q3q3 - _2q3 * m_ax + 4.0f * q0q0 * m_q1 - _2q0 * m_ay - _4q1 + _8q1 * q1q1 + _8q1 * q2q2 + _4q1 * m_az;
		s2 = 4.0f * q0q0 * m_q2 + _2q0 * m_ax + _4q2 * q3q3 - _2q3 * m_ay - _4q2 + _8q2 * q1q1 + _8q2 * q2q2 + _4q2 * m_az;
		s3 = 4.0f * q1q1 * m_q3 - _2q1 * m_ax + 4.0f * q2q2 * m_q3 - _2q2 * m_ay;
		recipNorm = invSqrt(s0 * s0 + s1 * s1 + s2 * s2 + s3 * s3); // normalize step magnitude
		s0 *= recipNorm;
		s1 *= recipNorm;
		s2 *= recipNorm;
		s3 *= recipNorm;

		// Apply feedback step
		qDot1 -= m_beta * s0;
		qDot2 -= m_beta * s1;
		qDot3 -= m_beta * s2;
		qDot4 -= m_beta * s3;
	}

	// Integrate rate of change of quaternion to yield quaternion
	m_q0 += qDot1 * m_dt_s;
	m_q1 += qDot2 * m_dt_s;
	m_q2 += qDot3 * m_dt_s;
	m_q3 += qDot4 * m_dt_s;

	// normalize quaternion
	recipNorm = invSqrt(m_q0 * m_q0 + m_q1 * m_q1 + m_q2 * m_q2 + m_q3 * m_q3);
	m_q0 *= recipNorm;
	m_q1 *= recipNorm;
	m_q2 *= recipNorm;
	m_q3 *= recipNorm;
}

// Fast inverse square-root
// See: http://en.wikipedia.org/wiki/Fast_inverse_square_root
float MADGWICK_AHRS::invSqrt(float x) {
	union {
		float f;
		uint32_t i;
	} conv = {x};	// member 'f' set to value of 'x'
		
	conv.i = 0x5f3759df - (conv.i >> 1);
	conv.f *= (1.5f - (0.5f * x * conv.f * conv.f));	// 1st iteration of the newton method
	conv.f *= (1.5f - (0.5f * x * conv.f * conv.f));	// 2nd iteration of the newton method (optional)
	
	return conv.f;
}
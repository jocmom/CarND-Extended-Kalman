#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/*
 * Constructor.
 */
FusionEKF::FusionEKF() {
  is_initialized_ = false;

  previous_timestamp_ = 0;

  // initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);
  Hj_ = MatrixXd(3, 4);

  //measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
        0, 0.0225;

  //measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
        0, 0.0009, 0,
        0, 0, 0.09;

  /**
  TODO:
    * Finish initializing the FusionEKF.
    * Set the process and measurement noises
    * The notation ν∼N(0,Q) defines the process noise as a
    * gaussian distribution with mean zero and covariance Q.
    *
    * The notation ω∼N(0,R) defines the measurement noise as a
    * gaussian distribution with mean zero and covariance R
    */
 	//measurement matrix
	H_laser_ << 1, 0, 0, 0,
    0, 1, 0, 0;
  //process noise
  noise_ax = 9;
  noise_ay = 9;
}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {


  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {
    /**
    TODO:
      * Initialize the state ekf_.x_ with the first measurement.
      * Create the covariance matrix.
      * Remember: you'll need to convert radar from polar to cartesian coordinates.
    */
    // first measurement
    cout << "EKF: " << endl;
    //state covariance matrix P
    Eigen::MatrixXd P = MatrixXd(4, 4);
    P << 1, 0, 0, 0,
      0, 1, 0, 0,
      0, 0, 1000, 0,
      0, 0, 0, 1000;
    //the initial transition matrix F_
    Eigen::MatrixXd F = MatrixXd(4, 4);
    F << 1, 0, 1, 0,
      0, 1, 0, 1,
      0, 0, 1, 0,
      0, 0, 0, 1;
    //the initial covariance matrix Q
    Eigen::MatrixXd Q = MatrixXd(4, 4);
    Q << 0, 0, 0, 0,
      0, 0, 0, 0,
      0, 0, 0, 0,
      0, 0, 0, 0;

    VectorXd x_in(4);
    x_in << 1, 1, 1, 1;

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
      cout << "Kalman Filter Initialization with RADAR values " << endl;
      //set the state with the initial location and velocity
      float rho = measurement_pack.raw_measurements_[0];
      float phi = measurement_pack.raw_measurements_[1];
      float rho_dot = measurement_pack.raw_measurements_[2];
      // Although radar gives velocity data in the form of the range rate
      // rho_dot, a radar measurement does not contain enough information
      // to determine the state variable velocities vx​ and v​y​​ 
      x_in << rho*cos(phi), rho*sin(phi), 0, 0; 
      ekf_.Init(x_in, P, F, Hj_, R_radar_, Q);
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */
      cout << "Kalman Filter Initialization with LIDAR values " << endl;
      //set the state with the initial location and zero velocity
      x_in << measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1], 0, 0;

      ekf_.Init(x_in, P, F, H_laser_, R_laser_, Q);
    }
    cout << x_in << endl;

    previous_timestamp_ = measurement_pack.timestamp_;
    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  /**
     * Update the state transition matrix F according to the new elapsed time.
      - Time is measured in seconds.
     * Update the process noise covariance matrix.
     * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
   */
	//compute the time elapsed between the current and previous measurements
  //dt - expressed in seconds
	float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;
	previous_timestamp_ = measurement_pack.timestamp_;

	//1. Modify the F matrix so that the time is integrated
	ekf_.F_(0, 2) = dt;
	ekf_.F_(1, 3) = dt;
	//2. Set the process covariance matrix Q
	float dt_2 = dt*dt;
	float dt_3 = dt_2*dt;
	float dt_4 = dt_3*dt;
	ekf_.Q_ << dt_4*noise_ax/4, 0, dt_3*noise_ax/2, 0,
    0, dt_4*noise_ay/4, 0, dt_3*noise_ay/2,
    dt_3*noise_ax/2, 0, dt_2*noise_ax, 0,
    0, dt_3*noise_ay/2, 0, dt_2*noise_ay;

	//3. Call the Kalman Filter predict() function
	//std::cout << "Predict" << endl;
  ekf_.Predict();

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  /**
   TODO:
     * Use the sensor type to perform the update step.
     * Update the state and covariance matrices.
   */

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar updates
    Tools tools;
    // for extended Kalman Filter the Jacobian matrix is required
    ekf_.H_ = tools.CalculateJacobian(ekf_.x_);
    ekf_.R_ = R_radar_;
    ekf_.UpdateEKF(measurement_pack.raw_measurements_);
  } else {
    // Laser updates
    ekf_.H_ = H_laser_; 
    ekf_.R_ = R_laser_;
    ekf_.Update(measurement_pack.raw_measurements_);
  }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}

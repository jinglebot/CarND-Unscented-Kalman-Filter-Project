#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 * This is scaffolding, do not modify
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 1.;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.5;
  
  //DO NOT MODIFY measurement noise values below these are provided by the sensor manufacturer.
  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;
  // std_radphi_ = 0.0175;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;
  // std_radrd_ = 0.1;

  //DO NOT MODIFY measurement noise values above these are provided by the sensor manufacturer.
  
  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */

  // initially set to false, set to true in first call of ProcessMeasurement
  is_initialized_ = false;

  // time when the state is true, in us
  time_us_ = 0;

  // Augmented state dimension
  n_x_ = 5;

  // Augmented state dimension
  n_aug_ = n_x_ + 2;

  //set weights
  lambda_ = 3 - n_aug_;
  weights_ = VectorXd(2 * n_aug_ + 1);
  weights_(0) = lambda_ / (lambda_ + n_aug_);
  for (int i = 1; i < 2 * n_aug_ + 1; i++) {
      weights_(i) = 0.5 / (lambda_ + n_aug_);
  }

  // NIS values
  NIS_radar = 0;
  NIS_laser = 0;

}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */

  cout << "time_us_: " << time_us_ << endl;

  // STATE DIMENSION INITIALIZATION
  if (!is_initialized_) {
    // state vector: [pos1 pos2 vel_abs yaw_angle yaw_rate] in SI units and rad
    x_.fill(0.0);  
    P_.fill(0.0);

    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {

      // Convert radar from polar to cartesian coordinates and initialize state.

      float rho = meas_package.raw_measurements_(0);
      float theta = meas_package.raw_measurements_(1);
      
      x_(0) = rho * cos (theta);
      x_(1) = rho * sin (theta);

      // initial covariance matrix
      P_ <<   std_radr_ * std_radr_,  0.,   0.,   0.,   0.,
              0.,   std_radphi_ * std_radphi_, 0.,   0.,   0., 
              0.,   0.,   std_radrd_ * std_radrd_,  0.,   0.,
              0.,   0.,   0.,   1.,   0.,
              0.,   0.,   0.,   0.,   1.;

    } else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {

      // Initialize state.
      x_(0) = meas_package.raw_measurements_(0);
      x_(1) = meas_package.raw_measurements_(1);

      // initial covariance matrix
      P_ <<   std_laspx_ * std_laspx_, 0.,   0.,   0.,   0.,
              0.,   std_laspy_ * std_laspy_, 0.,   0.,   0., 
              0.,   0.,   1.,   0.,   0.,
              0.,   0.,   0.,   1.,   0.,
              0.,   0.,   0.,   0.,   1.;

    }

    x_(2) = 4;
    x_(3) = 0.;
    x_(4) = 0.3;

    is_initialized_ = true;
    time_us_ = meas_package.timestamp_;
    return;
  }

  // PREDICTION
  double delta_t = (meas_package.timestamp_ - time_us_ ) / 1000000.0;
  time_us_ = meas_package.timestamp_;

  Prediction(delta_t);
  if ((use_radar_) && (meas_package.sensor_type_ == MeasurementPackage::RADAR)) UpdateRadar(meas_package);
  if ((use_laser_) && (meas_package.sensor_type_ == MeasurementPackage::LASER)) UpdateLidar(meas_package); 
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */

  // SIGMA POINTS CALCULATION

  // create sigma point matrix
  MatrixXd Xsig_ = MatrixXd(n_x_, 2 * n_x_ + 1);

  // calculate square root of P_
  MatrixXd A_ = P_.llt().matrixL();

  // Sigma point spreading parameter
  lambda_ = 3 - n_x_;

  //set sigma points as columns of matrix Xsig_
  Xsig_.col(0) = x_;
  for (int j = 0; j < n_x_; j++) {
    Xsig_.col(j + 1) = x_ + sqrt(lambda_ + n_x_) * A_.col(j);
    Xsig_.col(j + n_x_ + 1) = x_ - sqrt(lambda_ + n_x_) * A_.col(j);
  }

  // cout << "Xsig_ = " << endl << Xsig_ << endl;

  // AUGMENTED SIGMA POINTS CALCULATION

  //create augmented mean vector
  VectorXd x_aug_ = VectorXd(n_aug_);
  x_aug_.fill(0.0);
  x_aug_.head(n_x_) = x_;
  x_aug_(5) = 0.;
  x_aug_(6) = 0.;

  //create augmented state covariance
  MatrixXd P_aug_ = MatrixXd(n_aug_, n_aug_);
  P_aug_.fill(0.0);
  P_aug_.topLeftCorner(n_x_,n_x_) = P_;
  P_aug_(5,5) = std_a_ * std_a_;
  P_aug_(6,6) = std_yawdd_ * std_yawdd_;

  // calculate square root of P_aug_
  MatrixXd L_ = P_aug_.llt().matrixL();
  
  // Sigma point spreading parameter
  lambda_ = 3 - n_aug_;

  //create sigma point matrix
  MatrixXd Xsig_aug_ = MatrixXd(n_aug_, 2 * n_aug_ + 1);

  //create augmented sigma points
  Xsig_aug_.fill(0.0);
  Xsig_aug_.col(0)  = x_aug_;
  for (int i = 0; i< n_aug_; i++)
  {
    Xsig_aug_.col(i + 1) = x_aug_ + sqrt(lambda_ + n_aug_) * L_.col(i);
    Xsig_aug_.col(i + n_aug_ + 1) = x_aug_ - sqrt(lambda_ + n_aug_) * L_.col(i);
  }

  // cout << "Xsig_aug_ = " << endl << Xsig_aug_ << endl;
  
  // SIGMA POINTS PREDICTION

  Xsig_pred_ = MatrixXd(n_x_,2 * n_aug_ + 1);
  Xsig_pred_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {

    //extract values for better readability
    //double px = Xsig_aug_(0,i);
    //double py = Xsig_aug_(1,i);
    double v = Xsig_aug_(2,i);
    double psi = Xsig_aug_(3,i);
    double psi_dot = Xsig_aug_(4,i);
    double nu_ak = Xsig_aug_(5,i);
    double nu_psik = Xsig_aug_(6,i);
    
    // add noise
    VectorXd n_k1 = VectorXd(n_x_);
    n_k1 << 0.5 * delta_t * delta_t * cos(psi) * nu_ak,
            0.5 * delta_t * delta_t * sin(psi) * nu_ak,
            delta_t * nu_ak,
            0.5 * delta_t * delta_t * nu_psik,
            delta_t * nu_psik; 

    //avoid division by zero
    if (fabs(psi_dot) > 0.001) {
      VectorXd x_k1 = VectorXd(n_x_);
      x_k1 << (v / psi_dot) * (sin(psi + psi_dot * delta_t) - sin(psi)),
              (v / psi_dot) * (-cos(psi + psi_dot * delta_t) + cos(psi)),
              0,
              psi_dot * delta_t,
              0;

      Xsig_pred_.col(i) = x_ + x_k1 + n_k1;

    } else {
      VectorXd x_k1 = VectorXd(n_x_);
      x_k1 << v  * cos(psi) * delta_t,
              v * sin(psi) * delta_t,
              0,
              psi_dot * delta_t,
              0;

      Xsig_pred_.col(i) = x_ + x_k1 + n_k1;

    }
  }

  // cout << "Xsig_pred_ = " << endl << Xsig_pred_ << endl;

  // STATE MEAN AND COVARIANCE PREDICTION

  //predict state mean
  x_.fill(0.);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    x_ = x_ + weights_(i) * Xsig_pred_.col(i);
  }

  //predict state covariance matrix
  P_.fill(0.);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    VectorXd x_diff = Xsig_pred_.col(i) - x_;

    //angle normalization
    while (x_diff(3) > M_PI) x_diff(3) -= 2. * M_PI;
    while (x_diff(3) < -M_PI) x_diff(3) += 2. * M_PI;

    P_ = P_ + weights_(i) * x_diff * x_diff.transpose();
  }

  cout << "x_: " << x_ << "\n";
  cout << "P_: " << P_ << "\n";
} 


/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */
  //set radar measurement dimension: rho, phi, and rho_dot
  int n_z_ = 2;

  // Sigma point spreading parameter
  lambda_ = 3 - n_aug_;

  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z_, 2 * n_aug_ + 1);
  Zsig.fill(0.0);

  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z_);
  z_pred.fill(0.0);
  
  //measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z_,n_z_);
  S.fill(0.0);

  // SIGMA POINTS IN MEASUREMENT SPACE CALCULATION

  //transform sigma points into measurement space
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    double px = Xsig_pred_(0,i);
    double py = Xsig_pred_(1,i);
  
    VectorXd z = VectorXd(n_z_);
    z << px, py;
    Zsig.col(i) = z;
  }

  // calculate mean predicted measurement
  for (int m = 0; m < 2 * n_aug_ + 1; m++) {
    z_pred += Zsig.col(m) * weights_(m);
  }

  //calculate innovation covariance matrix S
  MatrixXd R = MatrixXd(n_z_,n_z_);
  R << std_laspx_ * std_laspx_,   0,
       0,                         std_laspy_ * std_laspy_;
  
  for (int n = 0; n < 2 * n_aug_ + 1; n++) {
    VectorXd z1 = VectorXd(n_z_);
    z1 = Zsig.col(n) - z_pred;
    S = S + weights_(n) * z1 * z1.transpose();
  }

  S = S + R;

  //create vector for incoming lidar measurement
  VectorXd z_in = VectorXd(n_z_);
  z_in = meas_package.raw_measurements_;

  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z_);
  Tc.fill(0.0);

  //calculate cross correlation matrix
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  

    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    
    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;

    // Angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }

  // Kalman gain K;
  MatrixXd K = Tc * S.inverse();

  // residual
  VectorXd z_diff = z_in - z_pred;

  // update state mean and covariance matrix
  x_ = x_ + K * z_diff;
  P_ = P_ - K * S * K.transpose();

  NIS_laser = z_in.transpose() * S.inverse() * z_in;

  cout << "x_: " << x_ << endl;
  cout << "P_: " << P_ << endl;
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */

  //set radar measurement dimension: rho, phi, and rho_dot
  int n_z_ = 3;

  // Sigma point spreading parameter
  lambda_ = 3 - n_aug_;

  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z_, 2 * n_aug_ + 1);
  Zsig.fill(0.0);

  // SIGMA POINTS IN MEASUREMENT SPACE CALCULATION

  //transform sigma points into measurement space
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    double px = Xsig_pred_(0,i);
    double py = Xsig_pred_(1,i);
    double v = Xsig_pred_(2,i);
    double psi = Xsig_pred_(3,i);
    //double psi_dot = Xsig_pred_(4,i);
    
    double rho = sqrt(px * px + py * py);
    double theta = atan2(py,px);
    double rho_dot;
    if (rho > 0.001)
          rho_dot = (px * cos(psi) * v + py * sin(psi) * v)/rho;
    else
          return;

    VectorXd z = VectorXd(n_z_);
    z << rho, theta, rho_dot;

    Zsig.col(i) = z;
  }

  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z_);
  z_pred.fill(0.0);
  
  // calculate mean predicted measurement
  for (int m = 0; m < 2 * n_aug_ + 1; m++) {
      z_pred += Zsig.col(m) * weights_(m);
  }

  // Angle normalization
  // while (z_pred(1) > M_PI) z_pred(1) -= 2. * M_PI;
  // while (z_pred(1) < -M_PI) z_pred(1) += 2. * M_PI;

  //calculate innovation covariance matrix S
  MatrixXd R = MatrixXd(n_z_,n_z_);
  R <<  std_radr_ * std_radr_,  0,                          0,
        0,                      std_radphi_ * std_radphi_,  0,
        0,                      0,                          std_radrd_ * std_radrd_;
  
  //measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z_,n_z_);
  S.fill(0.0);

  for (int n =0; n < 2 * n_aug_ + 1; n++) {
    //residual
    VectorXd z_diff = Zsig.col(n) - z_pred;

    // Angle normalization
    while (z_diff(1) > M_PI) z_diff(1) -= 2. * M_PI;
    while (z_diff(1) < -M_PI) z_diff(1) += 2. * M_PI;

    S = S + weights_(n) * z_diff * z_diff.transpose();
  }

  S = S + R;

  //cout << "z_pred: " << z_pred << endl;
  //cout << "S: " << S << endl;

  //create vector for incoming radar measurement
  VectorXd z_in = VectorXd(n_z_);
  z_in = meas_package.raw_measurements_;

  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z_);
  Tc.fill(0.0);

  //calculate cross correlation matrix
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  

    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    
    // Angle normalization
    while (z_diff(1) > M_PI) z_diff(1) -= 2. * M_PI;
    while (z_diff(1) < -M_PI) z_diff(1) += 2. * M_PI;

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;

    // Angle normalization
    while (x_diff(3) > M_PI) x_diff(3) -= 2. * M_PI;
    while (x_diff(3) < -M_PI) x_diff(3) += 2. * M_PI;

    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }

  // Kalman gain K;
  MatrixXd K = Tc * S.inverse();

  // residual
  VectorXd z_diff = z_in - z_pred;

  // Angle normalization
  while (z_diff(1) > M_PI) z_diff(1) -= 2. * M_PI;
  while (z_diff(1) < -M_PI) z_diff(1) += 2. * M_PI;

  // update state mean and covariance matrix
  x_ = x_ + K * z_diff;
  P_ = P_ - K * S * K.transpose();

  NIS_radar = z_in.transpose() * S.inverse() * z_in;
  
  cout << "x_: " << x_ << endl;
  cout << "P_: " << P_ << endl;
}

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
UKF::UKF(bool use_laser, bool use_radar) {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = use_laser;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = use_radar;

  // initial state vector
  x_ = VectorXd(n_x_);

  // initial covariance matrix
  P_ = MatrixXd(n_x_, n_x_);
  
  // Matrix to hold sigma points
  Xsig_pred_ = MatrixXd(n_x_, n_aug_size);
  
  // Vector for weights
  weights_ = VectorXd(n_aug_size);

  laser_R << std_laspx_ * std_laspx_, 0,
             0, std_laspy_ * std_laspy_;

  radar_R << std_radr_ * std_radr_, 0, 0,
             0, std_radphi_ * std_radphi_, 0,
             0, 0, std_radrd_ * std_radrd_;
  
}

UKF::~UKF() {}

bool UKF::isSupported(MeasurementPackage & meas_package) 
{
  if ( meas_package.sensor_type_ == MeasurementPackage::RADAR ) return use_radar_;
  if ( meas_package.sensor_type_ == MeasurementPackage::LASER ) return use_laser_;
  return false;
}

bool UKF::init(MeasurementPackage & meas_package) 
{
  if ( time_us_  < 0 ) {
    switch (meas_package.sensor_type_) {
      case MeasurementPackage::RADAR:
        radarInit(meas_package);
        time_us_ = meas_package.timestamp_;
        break;
      case MeasurementPackage::LASER:
        lidarInit(meas_package);
        time_us_ = meas_package.timestamp_;
        break;
    }
    return false;
  } else {
    // Calculate delta_t, store current time for future
    delta_t = (meas_package.timestamp_ - time_us_) / 1000000.0;
    time_us_ = meas_package.timestamp_;
    return true;
  }
}

void UKF::radarInit(MeasurementPackage & meas_package) 
{
      double rho = meas_package.raw_measurements_(0);
      double phi = meas_package.raw_measurements_(1);
      double rhodot = meas_package.raw_measurements_(2);
      // polar to cartesian - r * cos(angle) for x and r * sin(angle) for y**
      x_ << rho * cos(phi), rho * sin(phi), 4, rhodot * cos(phi), rhodot * sin(phi);
      //state covariance matrix
      P_ << std_radr_ * std_radr_, 0, 0, 0, 0,
            0, std_radr_ * std_radr_, 0, 0, 0,
            0, 0, 1, 0, 0,
            0, 0, 0, std_radphi_, 0,
            0, 0, 0, 0, std_radphi_;
}

void UKF::lidarInit(MeasurementPackage & meas_package)
{
      // ***** Last three values below can be tuned *****
      x_ << meas_package.raw_measurements_(0), meas_package.raw_measurements_(1), 4, 0.5, 0.0;
      //state covariance matrix
      P_ << std_laspx_*std_laspx_, 0, 0, 0, 0,
            0, std_laspy_*std_laspy_, 0, 0, 0,
            0, 0, 1, 0, 0,
            0, 0, 0, 1, 0,
            0, 0, 0, 0, 1;
}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage & meas_package) 
{
  if ( ! isSupported(meas_package) ) return;
  if ( init(meas_package) ) {
    // Predict
    Prediction();
    // Measurement updates
    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      UpdateRadar(meas_package);
    } else {
      UpdateLidar(meas_package);
    }
  }
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction() 
{
  //create augmented mean vector
  VectorXd x_aug_ = VectorXd(n_aug_);
  //create augmented state covariance
  MatrixXd P_aug_ = MatrixXd(n_aug_, n_aug_);
  //create sigma point matrix
  MatrixXd Xsig_aug_ = MatrixXd(n_aug_, n_aug_size);
  
  //create augmented mean state
  x_aug_.head(5) = x_;
  x_aug_(5) = 0;
  x_aug_(6) = 0;
  
  //create augmented covariance matrix
  MatrixXd Q = MatrixXd(2,2);
  Q << std_a_ * std_a_, 0,
        0, std_yawdd_ * std_yawdd_;
  P_aug_.fill(0.0);
  P_aug_.topLeftCorner(5, 5) = P_;
  P_aug_.bottomRightCorner(2, 2) = Q;
  
  //create square root matrix
  MatrixXd A_aug = P_aug_.llt().matrixL();
  
  //create augmented sigma points
  Xsig_aug_.col(0) = x_aug_;
  for (int i = 0; i < n_aug_; i++) {
    Xsig_aug_.col(i + 1) = x_aug_ + std::sqrt(lambda_ + n_aug_ ) * A_aug.col(i);
    Xsig_aug_.col(i + 1 + n_aug_) = x_aug_ - std::sqrt(lambda_ + n_aug_) * A_aug.col(i);
  }
  
  //predict sigma points
  //set vectors for each part added to x
  VectorXd vec1 = VectorXd(5);
  VectorXd vec2 = VectorXd(5);
  
  for (int i = 0; i < n_aug_size; ++i) {
    VectorXd calc_col = Xsig_aug_.col(i);
  //  double px = calc_col(0);
  //  double py = calc_col(1);
    double v = calc_col(2);
    double yaw = calc_col(3);
    double yawd = calc_col(4);
    double v_aug = calc_col(5);
    double v_yawdd = calc_col(6);
    
    //original
    VectorXd orig = calc_col.head(5);
    
    if (yawd > .001) {
      // If yaw dot is not zero
      vec1 << (v/yawd) * (sin(yaw + yawd * delta_t) - sin(yaw)),
              (v/yawd) * (-cos(yaw + yawd * delta_t) + cos(yaw)),
              0,
              yawd * delta_t,
              0;
    } else {
      // If yaw dot is zero - avoid division by zero
      vec1 << v * cos(yaw) * delta_t,
              v * sin(yaw) * delta_t,
              0,
              yawd * delta_t,
              0;
    }
    
    // This portion stays the same
    vec2 << .5 * delta_t * delta_t * cos(yaw) * v_aug,
            .5 * delta_t * delta_t * sin(yaw) * v_aug,
            delta_t * v_aug,
            .5 * delta_t * delta_t * v_yawdd,
            delta_t * v_yawdd;
    
    //write predicted sigma points into right column
    Xsig_pred_.col(i) << orig + vec1 + vec2;
  }
  
  //create vector for predicted state
  VectorXd x_pred = VectorXd(n_x_);
  //create covariance matrix for prediction
  MatrixXd P_pred = MatrixXd(n_x_, n_x_);
  
  x_pred.fill(0.0);
  P_pred.fill(0.0);
  
  for (int i = 0; i < n_aug_size; ++i) {
    //set weights
    weights_(i) = ( i == 0 ) ? lambda_ / (lambda_ + n_aug_) : .5 / (lambda_ + n_aug_);
    //predict state mean
    x_pred += weights_(i) * Xsig_pred_.col(i);
  }
  
  for (int i = 0; i < n_aug_size; ++i) {
    //predict state covariance matrix
    VectorXd x_diff = Xsig_pred_.col(i) - x_pred;
    //normalize angles
    normalizeAngles(3, x_diff);
    P_pred += weights_(i) * x_diff * x_diff.transpose();
  }
  
  x_ = x_pred;
  P_ = P_pred;
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage & meas_package) {  
  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(laser_dim, n_aug_size);
  //mean predicted measurement
  VectorXd z_pred = VectorXd(laser_dim);
  
  Zsig.fill(0.0);
  z_pred.fill(0.0);
  
  for (int i = 0; i < n_aug_size; ++i) {
    //transform sigma points into measurement space
    VectorXd state_vec = Xsig_pred_.col(i);
    Zsig.col(i) << state_vec(0),
                   state_vec(1);
    //calculate mean predicted measurement
    z_pred += weights_(i) * Zsig.col(i);
  }
  
  //measurement covariance matrix S
  MatrixXd S = calculateMeasurementCovarianceMatrix(Zsig, z_pred, laser_R, laser_dim, false);
  
  //create vector for incoming radar measurement
  VectorXd z = VectorXd(laser_dim);
  
  z << meas_package.raw_measurements_(0), //measurements px
       meas_package.raw_measurements_(1); //measurements py
  
  //create matrix for cross correlation Tc
  MatrixXd Tc = calculateCrossCorrelationMatrix(Zsig, z_pred, laser_dim, false);
  updateState(Tc, S, z, z_pred, & NIS_laser_, false);
  //cout<<"NIS laser is "<<NIS_laser_<<endl;
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage & meas_package) {  
  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(radar_dim, n_aug_size);
  
  //mean predicted measurement
  VectorXd z_pred = VectorXd(radar_dim);
  
  Zsig.fill(0.0);
  z_pred.fill(0.0);
  
  for (int i = 0; i < n_aug_size; ++i) {
    //transform sigma points into measurement space
    VectorXd state_vec = Xsig_pred_.col(i);
    double px = state_vec(0);
    double py = state_vec(1);
    double v = state_vec(2);
    double yaw = state_vec(3);
//    double yaw_d = state_vec(4);
    
    double rho = sqrt(px * px + py * py);
    double phi = atan2(py, px);
    double rho_d = (px * cos(yaw) * v + py * sin(yaw) *v) / rho;
    
    Zsig.col(i) << rho,
                   phi,
                   rho_d;
    
    //calculate mean predicted measurement
    z_pred += weights_(i) * Zsig.col(i);
  }
  
  //measurement covariance matrix S
  MatrixXd S = calculateMeasurementCovarianceMatrix(Zsig, z_pred, radar_R, radar_dim, true);
  
  //create example vector for incoming radar measurement
  VectorXd z = VectorXd(radar_dim);
  z << meas_package.raw_measurements_(0), //measurements rho
       meas_package.raw_measurements_(1), //measurements phi
       meas_package.raw_measurements_(2);  //measurements rhod
  
  //create matrix for cross correlation Tc
  MatrixXd Tc = calculateCrossCorrelationMatrix(Zsig, z_pred, radar_dim, true);
  updateState(Tc, S, z, z_pred, &NIS_radar_, true);
  //cout<<"NIS radar is "<<NIS_radar_<<endl;
}

MatrixXd UKF::calculateMeasurementCovarianceMatrix(MatrixXd & sig, VectorXd & pred, MatrixXd & R, int size, bool normalize)
{
  MatrixXd S = MatrixXd(size, size);
  S.fill(0.0);
  for (int i = 0; i < n_aug_size; i++) {
    VectorXd diff = sig.col(i) - pred;
    //normalize angles
    if (normalize) normalizeAngles(1, diff);
    S += weights_(i) * diff * diff.transpose();
  }
  // Add R to S
  S += R;
  return S;
}

MatrixXd UKF::calculateCrossCorrelationMatrix(MatrixXd & Zsig, VectorXd & z_pred, int size, bool normalize) 
{
  MatrixXd Tc = MatrixXd(n_x_, size);
  Tc.fill(0.0);
  for (int i = 0; i < n_aug_size; ++i) {
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    //normalize angles
    normalizeAngles(3, x_diff);
    VectorXd z_diff = Zsig.col(i) - z_pred;
    //normalize angles
    if (normalize) normalizeAngles(1, z_diff);
    Tc += weights_(i) * x_diff * z_diff.transpose();
  }
  return Tc;
}

void UKF::updateState(MatrixXd & Tc, MatrixXd & S, VectorXd & z, VectorXd & z_pred, double * nis, bool normalize) 
{
  // residual
  VectorXd z_diff = z - z_pred;
  //normalize angles
  if ( normalize ) normalizeAngles(1, z_diff);
  //calculate NIS
  *nis = z_diff.transpose() * S.inverse() * z_diff;

  //calculate Kalman gain K;
  MatrixXd K = Tc * S.inverse();
  //update state mean and covariance matrix
  x_ += K * z_diff;
  P_ -= K * S * K.transpose();
}

void UKF::normalizeAngles(int i, VectorXd & vect)
{
  if (vect(i) > M_PI) {
    vect(i) -= 2. * M_PI;
  } else if (vect(i) < -M_PI) {
    vect(i) += 2. * M_PI;
  }
}
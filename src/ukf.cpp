#include <cassert>
#include <iostream>
#include "ukf.h"
#include "Eigen/Dense"

using Eigen::MatrixXd;
using Eigen::VectorXd;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd::Zero(5);

  // initial covariance matrix
  // Let's take a diffuse prior so no worries about initialization
  P_ = 1e6*MatrixXd::Identity(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 30;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 30;
  
  /**
   * DO NOT MODIFY measurement noise values below.
   * These are provided by the sensor manufacturer.
   */

  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;
  
  /**
   * End DO NOT MODIFY section for measurement noise values 
   */
  
  /**
   * TODO: Complete the initialization. See ukf.h for other member properties.
   * Hint: one or more values initialized above might be wildly off...
   */
}

UKF::~UKF() {}

void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
    // Prediction first
    if(time_us_ >= 0){
        assert(meas_package.timestamp_ >= time_us_);
        double delta_t = 1.0e-6* (double)(meas_package.timestamp_ - time_us_);
        Prediction(delta_t);
    }
    time_us_ = meas_package.timestamp_;

    // And now measurement!
    if(meas_package.sensor_type_ == MeasurementPackage::RADAR)
        UpdateRadar(meas_package);

    if(meas_package.sensor_type_ == MeasurementPackage::LASER)
        UpdateLidar(meas_package);
}

void UKF::Prediction(double delta_t) {
  /**
   * TODO: Complete this function! Estimate the object's location. 
   * Modify the state vector, x_. Predict sigma points, the state, 
   * and the state covariance matrix.
   */
}

void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Use lidar data to update the belief 
   * about the object's position. Modify the state vector, x_, and 
   * covariance, P_.
   * You can also calculate the lidar NIS, if desired.
   */
}

void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Use radar data to update the belief 
   * about the object's position. Modify the state vector, x_, and 
   * covariance, P_.
   * You can also calculate the radar NIS, if desired.
   */
}


/**
 * Programming assignment functions:
 */


void GenerateSigmaPointsFunction(MatrixXd* Xsig_out, double lambda, VectorXd& x, MatrixXd& P ) {

    // set state dimension
    int n_x = x.rows();

    // calculate square root of P
    MatrixXd A = P.llt().matrixL();

    A*=sqrt(n_x + lambda);

    MatrixXd Xsig = MatrixXd(n_x, 2*n_x +1);
    Xsig.col(0) = x;
    for(int i=0;i<n_x; i++){
        Xsig.col(i+1) = x + A.col(i);
        Xsig.col(i+1+n_x)= x -A.col(i);
    }

    // print result
     std::cout << "Xsig = " << std::endl << Xsig << std::endl;
     Xsig_out = & Xsig;
}

template<class T>
void predictSingleVector(T& x, T& x_out, double dt){
    // let's give all elements of x_in names for readability
    double px = x[0];
    double py = x[1];
    double v = x[2];
    double psi = x[3];
    double psid = x[4];
    double nua = x[5];
    double nupsi = x[6];

    double  eps = 1e-6;

    // fill in x's
    if(fabs(psid) > eps){
        x_out[0] = v*(sin(psi + dt*psid)-sin(psi))/psid;
        x_out[1] = v*(-cos(psi + dt*psid)+cos(psi))/psid;
    }else{
        x_out[0] = v*cos(psi)*dt;
        x_out[1] = v*sin(psi)*dt;
    }
    x_out[0] += x[0] + 0.5*dt*dt*cos(psi)*nua;
    x_out[1] += x[1] + 0.5*dt*dt*sin(psi)*nua;

    x_out[2] = x[2] + dt*nua;
    x_out[3] = x[3] + psid*dt + 0.5*dt*dt*nupsi;
    x_out[4] = x[4] + dt*nupsi;

}

void UKF::GenerateSigmaPoints(MatrixXd* Xsig){
    GenerateSigmaPointsFunction(Xsig, 3-n_x_, x_, P_);
}

void UKF::AugmentedSigmaPoints(MatrixXd* Xsig_aug){

    // create augmented mean vector
    VectorXd x_aug = VectorXd::Zero(n_aug_);

    // create augmented state covariance
    MatrixXd P_aug = MatrixXd::Zero(n_aug_, n_aug_);

    /**
     * Student part begin
     */

    // create augmented mean state
    for(int i=0; i<n_x_; i++){
        x_aug(i) = x_(i);
    }

    // create augmented covariance matrix
    P_aug.block(0,0,5,5) = P_;
    P_aug(5,5) = std_a_*std_a_;
    P_aug(6,6) = std_yawdd_*std_yawdd_;
    // create square root matrix

    GenerateSigmaPointsFunction(Xsig_aug, 3-n_aug_, x_aug, P_aug);
}


void UKF::SigmaPointPrediction(MatrixXd* Xsig_out, MatrixXd& Xsig_aug) {
    // create matrix with predicted sigma points as columns
    MatrixXd Xsig_pred = MatrixXd(n_x_, 2 * n_aug_ + 1);

    double dt = 0.1; // time diff in sec

    // predict sigma points
    for(int i=0; i<Xsig_aug.cols(); i++){
        auto tmp1 = Xsig_aug.col(i);
        auto tmp2 = Xsig_pred.col(i);
        predictSingleVector(tmp1, tmp2, dt);
        //Xsig_pred.col(i) = tmp2;
    }

    // print result
    std::cout << "Xsig_aug = " << std::endl << Xsig_aug << std::endl;
    std::cout << "Xsig_pred = " << std::endl << Xsig_pred << std::endl;

    // write result
    *Xsig_out = Xsig_pred;
}

void UKF::PredictMeanAndCovariance(VectorXd* x_out, MatrixXd* P_out, MatrixXd& Xsig_pred){
    // set weights
    VectorXd weights = VectorXd::Constant(2*n_aug_+1, 0.5/(lambda_ + n_aug_));
    weights[0] = lambda_/(lambda_ + n_aug_);
    // predict state mean
    auto x = Xsig_pred * weights;
    // predict state covariance matrix
    MatrixXd deltas = MatrixXd(Xsig_pred);
    deltas.colwise() -= x;
    MatrixXd P = deltas*weights.asDiagonal()*deltas.transpose();
    /**
     * Student part end
     */

    // print result
    std::cout << "Predicted state" << std::endl;
    std::cout << x << std::endl;
    std::cout << "Predicted covariance matrix" << std::endl;
    std::cout << P << std::endl;

    // write result
    *x_out = x;
    *P_out = P;
}


template<class T>
void transformToRadar(T& x, T&z){
    double px = x[0];
    double py = x[1];
    double v = x[2];
    double psi = x[3];

    z[0] = sqrt(px*px + py*py);
    z[1] = atan(py/px);
    z[2] = (px*cos(psi) + py*sin(psi))*v/z[0];
}

void UKF::PredictRadarMeasurement(VectorXd* z_out, MatrixXd* S_out, MatrixXd& Xsig_pred) {
    // set measurement dimension, radar can measure r, phi, and r_dot
    int n_z = 3;

    // TODO: factor the weights out
    VectorXd weights = VectorXd(2*n_aug_+1);
    double weight_0 = lambda_/(lambda_+n_aug_);
    double weight = 0.5/(lambda_+n_aug_);
    weights(0) = weight_0;

    for (int i=1; i<2*n_aug_+1; ++i) {
        weights(i) = weight;
    }

    // radar measurement noise standard deviation radius in m
    //double std_radr = 0.3;

    // radar measurement noise standard deviation angle in rad
    //double std_radphi = 0.0175;

    // radar measurement noise standard deviation radius change in m/s
    //double std_radrd = 0.1;

    MatrixXd R = MatrixXd::Zero(3,3);
    R(0,0) = std_radr_*std_radr_;
    R(1,1) = std_radphi_*std_radphi_;
    R(2,2) = std_radrd_*std_radrd_;

    // create matrix for sigma points in measurement space
    MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);

    VectorXd z_pred = VectorXd(n_z);
    // transform sigma points into measurement space
    for(int i=0; i<Xsig_pred.cols(); i++){
        VectorXd tmp1 = Xsig_pred.col(i);
        transformToRadar(tmp1, z_pred);
        Zsig.col(i) = z_pred;
    }
    // calculate mean predicted measurement
    z_pred = Zsig*weights;
    // calculate innovation covariance matrix S
    MatrixXd deltas = MatrixXd(Zsig);
    deltas.colwise() -= z_pred;
    MatrixXd S = deltas*weights.asDiagonal()*deltas.transpose() +R;

    /**
     * Student part end
     */

    // print result
    std::cout << "z_pred: " << std::endl << z_pred << std::endl;
    std::cout << "S: " << std::endl << S << std::endl;

    // write result
    *z_out = z_pred;
    *S_out = S;
}
#include <cassert>
#include <iostream>
#include "ukf.h"
#include "Eigen/Dense"

using Eigen::MatrixXd;
using Eigen::VectorXd;


struct RadarSensorDetails: public SensorDetails{
public:
    RadarSensorDetails(){
        n_z = 3;

        // Radar measurement noise standard deviation radius in m
        double std_radr_ = 0.3;

        // Radar measurement noise standard deviation angle in rad
        double std_radphi_ = 0.03;

        // Radar measurement noise standard deviation radius change in m/s
        double std_radrd_ = 0.3;

        R = MatrixXd::Zero(3,3);
        R(0,0) = std_radr_*std_radr_;
        R(1,1) = std_radphi_*std_radphi_;
        R(2,2) = std_radrd_*std_radrd_;
    }

    void stateToMeasurement(VectorXd& x, VectorXd&z){
        double px = x[0];
        double py = x[1];
        double v = x[2];
        double psi = x[3];

        z[0] = sqrt(px*px + py*py);
        z[1] = atan(py/px);
        z[2] = (px*cos(psi) + py*sin(psi))*v/z[0];
    }

    ~RadarSensorDetails(){};

};

struct LidarSensorDetails: public SensorDetails{
public:
    LidarSensorDetails(){
        n_z = 2;

          // Laser measurement noise standard deviation position1 in m
  double std_laspx_=0.15;

  // Laser measurement noise standard deviation position2 in m
  double std_laspy_=0.15;


        R = MatrixXd::Zero(2,2);
        R(0,0) = std_laspx_*std_laspx_;
        R(1,1) = std_laspy_*std_laspy_;
    }

    void stateToMeasurement(VectorXd& x, VectorXd&z){
        double px = x[0];
        double py = x[1];

        z[0] = px;
        z[1] = py;
    }

    ~LidarSensorDetails(){};

};

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
  P_ = 10*MatrixXd::Identity(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 3;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 1;
  
  /**
   * DO NOT MODIFY measurement noise values below.
   * These are provided by the sensor manufacturer.
   */

// THESE WERE MOVED TO THE SensorDetails CLASSES ABOVE
//  // Laser measurement noise standard deviation position1 in m
//  std_laspx_ = 0.15;
//
//  // Laser measurement noise standard deviation position2 in m
//  std_laspy_ = 0.15;

//  // Radar measurement noise standard deviation radius in m
//  std_radr_ = 0.3;
//
//  // Radar measurement noise standard deviation angle in rad
//  std_radphi_ = 0.03;
//
//  // Radar measurement noise standard deviation radius change in m/s
//  std_radrd_ = 0.3;
  
  /**
   * End DO NOT MODIFY section for measurement noise values 
   */
  
  /**
   * TODO: Complete the initialization. See ukf.h for other member properties.
   * Hint: one or more values initialized above might be wildly off...
   */
   n_x_ = 5;
   n_aug_=7;
    lambda_ = 3-n_aug_;

    weights_ = VectorXd::Constant(2*n_aug_+1, 0.5/(lambda_ + n_aug_));
    weights_[0] = lambda_/(lambda_ + n_aug_);

    Xsig_aug_ = MatrixXd(n_aug_, 2*n_aug_ + 1);
    Xsig_pred_ = MatrixXd(n_x_, 2*n_aug_ + 1);
}

UKF::~UKF() {}

void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
    // Prediction first
    double delta_t;
    if(time_us_ >= 0){
        assert(meas_package.timestamp_ >= time_us_);
        delta_t = 1.0e-6* (double)(meas_package.timestamp_ - time_us_);

    }else{
        delta_t = 0.0;
    }
    time_us_ = meas_package.timestamp_;

    Prediction(delta_t);

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
  // calculate the augmented sigma points
   AugmentedSigmaPoints();
   //Predict the sigma points
   SigmaPointPrediction(delta_t);
   PredictMeanAndCovariance();

}

void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Use lidar data to update the belief 
   * about the object's position. Modify the state vector, x_, and 
   * covariance, P_.
   * You can also calculate the lidar NIS, if desired.
   */
    auto sd = LidarSensorDetails();
    ProcessSensor(meas_package.raw_measurements_,// actual measurement
                  sd
    );
}

void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Use radar data to update the belief 
   * about the object's position. Modify the state vector, x_, and 
   * covariance, P_.
   * You can also calculate the radar NIS, if desired.
   */
   auto sd = RadarSensorDetails();
    ProcessSensor(meas_package.raw_measurements_,// actual measurement
                  sd
    );
}


/**
 * Programming assignment functions:
 */


void GenerateSigmaPointsFunction(MatrixXd& Xsig_out, double lambda, VectorXd& x, MatrixXd& P ) {

    // set state dimension
    int n_x = x.rows();

    // calculate square root of P
    MatrixXd A = P.llt().matrixL();

    A*=sqrt(n_x + lambda);

    Xsig_out.col(0) = x;
    for(int i=0;i<n_x; i++){
        Xsig_out.col(i+1) = x + A.col(i);
        Xsig_out.col(i+1+n_x)= x -A.col(i);
    }

//    // print result
//     std::cout << "Xsig = " << std::endl << Xsig << std::endl;
//     Xsig_out = & Xsig;
}

void UKF::AugmentedSigmaPoints(){
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

    GenerateSigmaPointsFunction(Xsig_aug_, 3-n_aug_, x_aug, P_aug);
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

void UKF::SigmaPointPrediction(double dt) {
    // create matrix with predicted sigma points as columns
    //MatrixXd Xsig_pred = MatrixXd(n_x_, 2 * n_aug_ + 1);
    // predict sigma points
    for(int i=0; i<Xsig_aug_.cols(); i++){
        auto tmp1 = Xsig_aug_.col(i);
        auto tmp2 = Xsig_pred_.col(i);
        predictSingleVector(tmp1, tmp2, dt);
        Xsig_pred_.col(i) = tmp2;
    }

//    // print result
//    std::cout << "Xsig_aug = " << std::endl << Xsig_aug << std::endl;
//    std::cout << "Xsig_pred = " << std::endl << Xsig_pred << std::endl;
//
//    // write result
//    *Xsig_out = Xsig_pred;
}

void UKF::PredictMeanAndCovariance(){
    // set weights

    // predict state mean
    x_ = Xsig_pred_ * weights_;
    // predict state covariance matrix
    X_deltas_ = MatrixXd(Xsig_pred_);
    X_deltas_.colwise() -= x_;
    P_ = X_deltas_*weights_.asDiagonal()*X_deltas_.transpose();

//    // print result
//    std::cout << "Predicted state" << std::endl;
//    std::cout << x_ << std::endl;
//    std::cout << "Predicted covariance matrix" << std::endl;
//    std::cout << P_ << std::endl;
}


void UKF::ProcessSensor(VectorXd& z,
        SensorDetails& sd// actual measurement
        ) {
    // create matrix for sigma points in measurement space
    MatrixXd Zsig = MatrixXd(sd.n_z, 2 * n_aug_ + 1);
    VectorXd z_pred = VectorXd(sd.n_z);

    // transform sigma points into measurement space
    for(int i=0; i<Xsig_pred_.cols(); i++){
        VectorXd tmp1 = Xsig_pred_.col(i);
        sd.stateToMeasurement(tmp1, z_pred);
        Zsig.col(i) = z_pred;
    }
    // calculate mean predicted measurement
    z_pred = Zsig*weights_;
    // calculate innovation covariance matrix S
    MatrixXd Z_deltas = MatrixXd(Zsig);
    Z_deltas.colwise() -= z_pred;
    MatrixXd S = Z_deltas*weights_.asDiagonal()*Z_deltas.transpose() + sd.R;

//    // print result
//    std::cout << "z_pred: " << std::endl << z_pred << std::endl;
//    std::cout << "S: " << std::endl << S << std::endl;

    // HERE THE MERGE FROM THE OTHER FUNCTION
    MatrixXd Tc = X_deltas_*weights_.asDiagonal()*Z_deltas.transpose();
    MatrixXd K = Tc*S.inverse();

    x_ += K *(z - z_pred);
    P_ -= K*S*K.transpose();
}


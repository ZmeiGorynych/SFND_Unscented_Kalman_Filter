#ifndef UKF_H
#define UKF_H

#include "Eigen/Dense"
#include "measurement_package.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

struct SensorDetails {
public:
    SensorDetails() {};

    virtual ~SensorDetails() {};

    int n_z; // dimensionality of sensor output
    Eigen::MatrixXd R; // measurement noise
    // converts a single state vector to the corresponding observation
    virtual void stateToMeasurement(VectorXd &x, VectorXd &z) = 0;
};


class UKF {
public:
    /**
     * Constructor
     */
    UKF();

    /**
     * Destructor
     */
    virtual ~UKF();

    /**
     * ProcessMeasurement
     * @param meas_package The latest measurement data of either radar or laser
     */
    void ProcessMeasurement(MeasurementPackage meas_package);

    /**
     * Prediction Predicts sigma points, the state, and the state covariance
     * matrix
     * @param delta_t Time between k and k+1 in s
     */
    void Prediction(double delta_t);

    /**
     * Updates the state and the state covariance matrix using a laser measurement
     * @param meas_package The measurement at k+1
     */
    void UpdateLidar(MeasurementPackage meas_package);

    /**
     * Updates the state and the state covariance matrix using a radar measurement
     * @param meas_package The measurement at k+1
     */
    void UpdateRadar(MeasurementPackage meas_package);

    //void GenerateSigmaPoints(MatrixXd* Xsig);
    void AugmentedSigmaPoints();

    void SigmaPointPrediction(double dt);

    void PredictMeanAndCovariance();

    void ProcessSensor(VectorXd &z,// actual measurement
                       SensorDetails &sd,
                       bool error_mod_2pi
    );
    // initially set to false, set to true in first call of ProcessMeasurement
//  bool is_initialized_;

    // if this is false, laser measurements will be ignored (except for init)
    bool use_laser_;

    // if this is false, radar measurements will be ignored (except for init)
    bool use_radar_;

    // state vector: [pos1 pos2 vel_abs yaw_angle yaw_rate] in SI units and rad
    Eigen::VectorXd x_;

    // state covariance matrix
    Eigen::MatrixXd P_;

    // predicted sigma points matrix
    Eigen::MatrixXd Xsig_aug_, Xsig_pred_, X_deltas_;

    // time when the state is true, in us
    long long time_us_ = -1;

    // Process noise standard deviation longitudinal acceleration in m/s^2
    double std_a_;

    // Process noise standard deviation yaw acceleration in rad/s^2
    double std_yawdd_;

    // Weights of sigma points
    Eigen::VectorXd weights_;

    // State dimension
    int n_x_;

    // Augmented state dimension
    int n_aug_;

    // Sigma point spreading parameter
    double lambda_;
};


#endif  // UKF_H
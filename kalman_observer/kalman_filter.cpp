// Kalman filter class implementation
// @author Giada Simionato, 1822614

#include <iostream>
#include <stdexcept>
#include <math.h>
#include "Eigen/Dense"
#include "kalman_filter.hpp"
using namespace std;


// Constructor implementation
KalmanFilter::KalmanFilter(double dt, const Eigen::MatrixXd& A, const Eigen::MatrixXd& B, const Eigen::MatrixXd& C, const Eigen::MatrixXd& R, double sigma_jerk, double sigma_ddfext){

    Eigen::Matrix2d cov_input;
    cov_input << sigma_jerk, 0, sigma_ddfext, 0;
    Q = B*cov_input*B.transpose();                  // compute covariance process noise
    this->dt = dt;
    this->A = A;
    this->C = C;
    this->R = R;
    m = C.rows();
    n = A.rows();
    x_hat.resize(n);

}


// KF initialization no guessing
void KalmanFilter::init(){

    int f = pow(10, 2);
    x_hat.setZero();
    Eigen::MatrixXd p(n,n);
    P = p;
    P = f*P.setIdentity();
    t = 0;
}

// KF initialization with guessing
void KalmanFilter::init(const Eigen::VectorXd& x0, const Eigen::MatrixXd& P0){

    x_hat = x0;
    P = P0;
    t = 0;
}

// Update function implementation (for Z observer)
void KalmanFilter::update(const Eigen::VectorXd& y){

    // Prediction step
    x_hat = A*x_hat;
    P = A*P*A.transpose() + Q;

    // Correction step
    K = P*C.transpose()*(C*P*C.transpose()+R).inverse();
    NI = y - C*x_hat;
    x_hat = x_hat + K*NI;
    P = P - K*C*P;

    t += dt;    // Get next instant
}

// @overloaded update function implementation (for X, Y observers)
void KalmanFilter::update(const Eigen::MatrixXd C, const Eigen::VectorXd& y){
    this->C = C;
    update(y);
}

// Returns current state estimate
Eigen::VectorXd KalmanFilter::getState(){
    return x_hat;
}

// Returns current timestep
double KalmanFilter::getTimestep(){
    return t;
}

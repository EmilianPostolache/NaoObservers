// Kalman filter class implementation
// From paper "State Estimation for Force-Controlled Humanoid Balance using Simple Models in the Presence of Modeling Error"
// @author Riccardo Ratini, 1656801

#include <iostream>
#include <stdexcept>
#include <math.h>
#include "Eigen/Dense"
#include "ef_estimator.hpp"
using namespace std;


// Constructor implementation
ExtForceKalmanFilter::ExtForceKalmanFilter(double dt, const Eigen::MatrixXd& A, const Eigen::MatrixXd& B, const Eigen::MatrixXd& C, const Eigen::MatrixXd& Q, const Eigen::MatrixXd& R){

    this->dt = dt;
    this->A = A;
    this->B = B;
    this->C = C;
    this->Q = Q;
    this->R = R;
    m = C.rows();
    n = A.rows();
    x_hat.resize(n);

}


// KF initialization no guessing
void ExtForceKalmanFilter::init(){

    int f = pow(10, 2);
    x_hat.setZero();
    P.resize(n,n);
    P = f*P.setIdentity();
    t = 0;
}

// KF initialization with guessing
void ExtForceKalmanFilter::init(const Eigen::VectorXd& x0, const Eigen::MatrixXd& P0){

    x_hat = x0;
    P = P0;
    t = 0;
}

// Update function implementation
void ExtForceKalmanFilter::update(const Eigen::VectorXd& u, const Eigen::VectorXd& y){

    // Prediction step
    x_hat = A*x_hat + B*u;
    P = A*P*A.transpose() + Q;

    // Correction step
    G = P*C.transpose()*(C*P*C.transpose()+R).inverse();
    NI = y - C*x_hat;
    x_hat = x_hat + G*NI;
    P = P - G*C*P;

    t += dt;    // Get next instant
}


// Returns current state estimate
Eigen::VectorXd ExtForceKalmanFilter::getState(){
    return x_hat;
}

// Returns current timestep
double ExtForceKalmanFilter::getTimestep(){
    return t;
}

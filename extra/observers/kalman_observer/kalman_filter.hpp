// Kalman filter class header file
// @author Giada Simionato, 1822614

//#include <Eigen/Dense>
#include "Eigen/Dense"

class KalmanFilter {

public:

    /*
    Class constructor

    @param dt: timestep
    @param A: system dynamics matrix
    @param B: input matrix
    @param C: output matrix
    @param R: process noise covariance
    @param P: estimate error covariance
    @param sigma_jerk: covariance of CoM jerk
    @param sigma_ddfext: covariance of second derivative of external force
    */

    KalmanFilter(double dt,
                 const Eigen::MatrixXd& A,
                 const Eigen::MatrixXd& B,
                 const Eigen::MatrixXd& C,
                 const Eigen::MatrixXd& R,
                 double sigma_jerk,
                 double sigma_ddfext);


    /*
    Initialization of the Kalman Filter
    */

    void init();


    /*
    Initialization of the Kalman Filter with initial guessing

    @param x0: initial estimate of state system
    @param P0: initial covariance estimate
    */

    void init(const Eigen::VectorXd& x0, const Eigen::MatrixXd& P0);


    /*
    Update step (contains both prediction and correction) - for Z observer

    @param y: vector of measurements
    */
    
    void update(const Eigen::VectorXd& y);

    /*
    Update step (contains both prediction and correction) - for X, Y observers

    @param y: vector of measurements
    @param C: dynamic matrix to overwrite
    */

    void update(const Eigen::MatrixXd C, const Eigen::VectorXd& y);
    

    /*
    Returns the current estimate of the system state
    */

    Eigen::VectorXd getState();


    /*
    Returns the current timestep
    */

    double getTimestep();


    private:

        int m, n; // System dimensions
        Eigen::MatrixXd A, C, Q, R, P, K, NI; //Kalman Filter matrices
        Eigen::VectorXd x_hat; // state system estimate
        double t, dt;  // time elapsed, timestep

};
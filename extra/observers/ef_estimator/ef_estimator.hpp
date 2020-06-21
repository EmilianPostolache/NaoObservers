// Kalman filter class header file
// From paper "State Estimation for Force-Controlled Humanoid Balance using Simple Models in the Presence of Modeling Error"
// @author Riccardo Ratini, 1656801

//#include <Eigen/Dense>
#include "Eigen/Dense"

class ExtForceKalmanFilter {

public:

    /*
    Class constructor

    @param dt: timestep
    @param A: system dynamics matrix
    @param B: input matrix
    @param C: output matrix
    @param Q: process noice covariance 
    @param R: output noise covariance
    @param P: estimate error covariance
    */
    ExtForceKalmanFilter(double dt,
            const Eigen::MatrixXd& A,
            const Eigen::MatrixXd& B,
            const Eigen::MatrixXd& C,
            const Eigen::MatrixXd& Q,
            const Eigen::MatrixXd& R);


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
    Update step (contains both prediction and correction)

    @param x_ vectpr of inputs
    @param y: vector of measurements
    */
    void update(const Eigen::VectorXd& u, const Eigen::VectorXd& y);

    

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
        Eigen::MatrixXd A, B, C, Q, R, P, G, NI; //Kalman Filter matrices
        Eigen::VectorXd x_hat; // state system estimate
        double t, dt;  // time elapsed, timestep

};
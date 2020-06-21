// .cpp file for testing the ExtForceKalmanFilter

#include <iostream>
#include <fstream>
using namespace std;
#include <math.h>
#include "Eigen/Dense"
#include "Eigen/Eigen"
#include "ef_estimator.hpp"


// Function that retrieves the row with index r in the matrix m (no other shortcuts in the Eigen documentation worked)
Eigen::VectorXd getRow(int r, Eigen::MatrixXd m){

      Eigen::VectorXd row(m.cols());
      for (int i=0; i<m.cols(); i++){
            row(i) = m(r,i);
      }
      return row;
}


int main(){

    // Initialize values experiment
    
    int n = 4;          // size vector states
    int m = 2;          // size vector measurements
    int u = 1;          // size input vector
    int ny = 30;        // number measurements ----------------------------- CHANGE HERE NUMBER SAMPLES -----------------------------------

    double freq = 60;   // [Hz] sampling rate
    double dt = 1/freq; // [s] timestep

    double g = 9.81;    // [m/s^2] gravity acceleration
    double ni = sqrt(g/0.33);

    double pos_pro_noise = exp(-8), vel_pro_noise = exp(-4), output_noise = exp(-5);

    // Initialize matrices for observers

    Eigen::MatrixXd Q(n, n);
    Eigen::MatrixXd R(m, m);

    Eigen::MatrixXd A(n, n);
    Eigen::MatrixXd B(n, u);
    Eigen::MatrixXd C(m, n);

    // Define matrices for observers
    Q <<  pos_pro_noise, 0, 0, 0,
          0, vel_pro_noise, 0, 0,
          0, 0, vel_pro_noise, 0,
          0, 0, 0, vel_pro_noise;

    R << output_noise, 0,
         0, output_noise;

    A << 1,           dt,             0, 0,
         pow(ni,2)*dt, 1, -pow(ni,2)*dt, dt, 
         0,            0,             1, 0, 
         0,            0,             0, 1;
        
    B << 0, 0, 1, 0;

    C << 1, 0, 0, 0, 
         0, 0, 1, 0;

    // Initialization Observers and vector of measurements
    ExtForceKalmanFilter observer(dt, A, B, C, Q, R);  
    observer.init();

    Eigen::VectorXd u(ny);
    Eigen::VectorXd Xc(ny);
    Eigen::VectorXd Xz(ny);

    u << 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 
         0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 
         0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1;

    Xc << 0.09, 0.214, 0.32, 0.382, 0.512, 0.584, 0.72, 0.818, 0.89, 1.02, 
          1.084, 1.19, 1.308, 1.404, 1.486, 1.592, 1.69, 1.808, 1.88, 2.004, 
          2.118, 2.192, 2.288, 2.408, 2.498, 2.6, 2.716, 2.788, 2.884, 3.012;

    Xz << 0.098, 0.192, 0.28, 0.386, 0.51, 0.582, 0.708, 0.818, 0.892, 0.994, 
          1.092, 1.196, 1.31, 1.398, 1.52, 1.586, 1.702, 1.786, 1.898, 2.01, 
          2.104, 2.2, 2.306, 2.398, 2.516, 2.58, 2.688, 2.802, 2.91, 2.982;


     for (size_t i = 0; i < 30; i++)
     {
          Eigen::VectorXd u_i(1);
          Eigen::VectorXd y(2);
          u_i << u[i];
          y << Xc[i], Xz[i];

          std::cout << "Timestamp: "  << i << std::endl;
          std::cout << "Input u : "  << u_i << std::endl;
          std::cout << "Measures : "  << y << std::endl;
          observer.update(u_i, y);
          std::cout << "Updated State: " << observer.state() << std::endl;
          std::cout << "=====================================" << std::endl;
     }
    
    return 0;
}

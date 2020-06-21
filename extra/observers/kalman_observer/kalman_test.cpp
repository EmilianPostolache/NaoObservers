// .cpp file for testing the Kalman Filter

#include <iostream>
#include <fstream>
using namespace std;
#include <math.h>
#include "Eigen/Dense"
#include "Eigen/Eigen"
#include "kalman_filter.hpp"


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
      
      int n = 5;          // size vector states
      int m = 3;          // size vector measurements
      int u = 2;          // size input vector
      int ny = 30;        // number measurements ----------------------------- CHANGE HERE NUMBER SAMPLES -----------------------------------

      double freq = 60;   // [Hz] sampling rate
      double dt = 1/freq; // [s] timestep

      double g = 9.81;    // [m/s^2] gravity acceleration
      double Mc = 5.19;   // [kg] mass of the robot

      double sigmaQ_jerk = pow(10.0, 3.0);
      double sigmaQ_ddfext = pow(10.0, 3.0);

      double zc_hat, ddzc_hat, fz_hat, grf_hat;
      double fx_hat, fy_hat, ts;

      // Initialize matrices for observers

      Eigen::MatrixXd Rz(m, m);
      Eigen::MatrixXd Rx(m,m);
      Eigen::MatrixXd Ry(m,m);

      Eigen::MatrixXd A(n, n);
      Eigen::MatrixXd B(n, u);
      Eigen::MatrixXd Cz(m, n);
      Eigen::MatrixXd Cx(m, n);
      Eigen::MatrixXd Cy(m, n);

      // Define matrices for observers

      Rz << 0.01, 0, 0, 0, 1, 0, 0, 0, 1;
      Rx << 0.01, 0, 0, 0, 1, 0, 0, 0, 0.01;
      Ry = Rx;

      A << 1, dt, pow(dt, 2.0)/2, 0, 0,
            0, 1, dt, 0, 0, 
            0, 0, 1, 0, 0, 
            0, 0, 0, 1, dt, 
            0, 0, 0, 0, 1;
            
      B << pow(dt, 3.0)/6, 0, pow(dt, 2.0)/2, 0, dt, 0, 0, pow(dt, 2.0)/2, 0, dt;

      Cz << 1, 0, 0, 0, 0, 
            0, 0, 1, 0, 0, 
            0, 0, -Mc, 1, 0;

      // Initialization Observers and vector of measurements

      KalmanFilter z_observer(dt, A, B, Cz, Rz, sigmaQ_jerk, sigmaQ_ddfext);  // z-axis Observer
      z_observer.init();

      KalmanFilter x_observer(dt, A, B, Cz, Rx, sigmaQ_jerk, sigmaQ_ddfext);  // x-axis Observer
      x_observer.init();

      KalmanFilter y_observer(dt, A, B, Cz, Ry, sigmaQ_jerk, sigmaQ_ddfext);  // z-axis Observer
      y_observer.init();

      // ------------------------------------------<INSERT IN Y_z, Y_x, Y_y DART VALUES as described in rows 75%77> --------------------------------------

      Eigen::MatrixXd Y_z(ny, m);  // measurements matrix for z-observer. ROWS: measurements for each timestep, COLS: z_CoM, ddz_CoM, tot_FSR+Mc*g
      Eigen::MatrixXd Y_x(ny, m);  // measurements matrix for x-observer. ROWS: measurements for each timestep, COLS: x_CoM, ddx_CoM, x_ZMP
      Eigen::MatrixXd Y_y(ny, m);  // measurements matrix for y-observer. ROWS: measurements for each timestep, COLS: y_CoM, ddy_CoM, y_ZMP
      Eigen::VectorXd y_z(m);     // measurements for z-observer
      Eigen::VectorXd y_x(m);     // measurements for x-observer
      Eigen::VectorXd y_y(m);     // measurements for y-observer

      // DUMMY VALUES AS RICCARDO

      Y_z << 0.09, 0.214, 0.32, 0.382, 0.512, 0.584, 0.72, 0.818, 0.89, 1.02, 
            1.084, 1.19, 1.308, 1.404, 1.486, 1.592, 1.69, 1.808, 1.88, 2.004, 
            2.118, 2.192, 2.288, 2.408, 2.498, 2.6, 2.716, 2.788, 2.884, 3.012,
            0.09, 0.214, 0.32, 0.382, 0.512, 0.584, 0.72, 0.818, 0.89, 1.02, 
            1.084, 1.19, 1.308, 1.404, 1.486, 1.592, 1.69, 1.808, 1.88, 2.004, 
            2.118, 2.192, 2.288, 2.408, 2.498, 2.6, 2.716, 2.788, 2.884, 3.012,
            0.09, 0.214, 0.32, 0.382, 0.512, 0.584, 0.72, 0.818, 0.89, 1.02, 
            1.084, 1.19, 1.308, 1.404, 1.486, 1.592, 1.69, 1.808, 1.88, 2.004, 
            2.118, 2.192, 2.288, 2.408, 2.498, 2.6, 2.716, 2.788, 2.884, 3.012;

      Y_x << 0.09, 0.214, 0.32, 0.382, 0.512, 0.584, 0.72, 0.818, 0.89, 1.02, 
            1.084, 1.19, 1.308, 1.404, 1.486, 1.592, 1.69, 1.808, 1.88, 2.004, 
            2.118, 2.192, 2.288, 2.408, 2.498, 2.6, 2.716, 2.788, 2.884, 3.012,
            0.09, 0.214, 0.32, 0.382, 0.512, 0.584, 0.72, 0.818, 0.89, 1.02, 
            1.084, 1.19, 1.308, 1.404, 1.486, 1.592, 1.69, 1.808, 1.88, 2.004, 
            2.118, 2.192, 2.288, 2.408, 2.498, 2.6, 2.716, 2.788, 2.884, 3.012,
            0.09, 0.214, 0.32, 0.382, 0.512, 0.584, 0.72, 0.818, 0.89, 1.02, 
            1.084, 1.19, 1.308, 1.404, 1.486, 1.592, 1.69, 1.808, 1.88, 2.004, 
            2.118, 2.192, 2.288, 2.408, 2.498, 2.6, 2.716, 2.788, 2.884, 3.012;

      Y_y << 0.09, 0.214, 0.32, 0.382, 0.512, 0.584, 0.72, 0.818, 0.89, 1.02, 
            1.084, 1.19, 1.308, 1.404, 1.486, 1.592, 1.69, 1.808, 1.88, 2.004, 
            2.118, 2.192, 2.288, 2.408, 2.498, 2.6, 2.716, 2.788, 2.884, 3.012,
            0.09, 0.214, 0.32, 0.382, 0.512, 0.584, 0.72, 0.818, 0.89, 1.02, 
            1.084, 1.19, 1.308, 1.404, 1.486, 1.592, 1.69, 1.808, 1.88, 2.004, 
            2.118, 2.192, 2.288, 2.408, 2.498, 2.6, 2.716, 2.788, 2.884, 3.012,
            0.09, 0.214, 0.32, 0.382, 0.512, 0.584, 0.72, 0.818, 0.89, 1.02, 
            1.084, 1.19, 1.308, 1.404, 1.486, 1.592, 1.69, 1.808, 1.88, 2.004, 
            2.118, 2.192, 2.288, 2.408, 2.498, 2.6, 2.716, 2.788, 2.884, 3.012;

      // For each measurement

      ofstream ext_file("Force_estimates.txt");

      for (int r=0; r<ny; r++){

          y_z = getRow(r, Y_z);
          y_x = getRow(r, Y_x);
          y_y = getRow(r, Y_y);
    
          z_observer.update(y_z);

          zc_hat = z_observer.getState()(0);   // estimate of z CoM
          ddzc_hat = z_observer.getState()(2);   // estimate of acceleration z CoM
          fz_hat = z_observer.getState()(3);   // estimate of external force z axis
          grf_hat = -Mc*g -Mc*ddzc_hat + fz_hat;

          Cx << 1, 0, 0, 0, 0,
                0, 0, 1, 0, 0,
                1, 0, Mc*zc_hat/grf_hat, -zc_hat/grf_hat, 0;

          Cy = Cx;

          x_observer.update(Cx, y_x);
          fx_hat = x_observer.getState()(3);   // estimate of external force x axis
          x_observer.update(Cy, y_y);
          fy_hat = y_observer.getState()(3);   // estimate of external force y axis
          ts = z_observer.getTimestep();
          ext_file << ts << "," << fx_hat << "," << fy_hat << "," << fz_hat << endl;
    }
    ext_file.close();
    
    return 0;
}

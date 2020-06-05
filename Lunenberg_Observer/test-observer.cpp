#include <iostream>

#include <vector>
#include "Eigen/Dense"

#include "observer.hpp"

int main(int argc, char* argv[]){
    int n = 5;
    int m = 2;

    double dt = 1.0/30;

    Eigen::MatrixXd A(n, n);
    Eigen::MatrixXd B(n, 1);
    Eigen::MatrixXd C(m, n);
    Eigen::MatrixXd G(n, m);

    double ni = sqrt(9.81/0.33);

    A << 0, 1, 0, 0, 0, 
         pow(ni,2), 0, -pow(ni,2), 1, 0,
         0, 0, 0, 0, 0,
         0, 0, 0, 0, 1,
         0, 0, 0, 0, 0;

    B << 0, 0, 1, 0, 0;

    C << 1, 0, 0, 0, 0,
         0, 0, 1, 0, 0;

     G << 0, -1.36,
         29, -29.74,
         -0.05, 0,
         0, 0.68,
         0.07, 0.01;

    /*G << 11.0276, -2.6295,
         70.9776, -49.8010,
         -0.0214, 3.9724,
         61.5672, -42.7015,
         30.3446, -25.2731;
*/
    std::cout   <<  "A: \n" << A << std::endl;
    std::cout   <<  "B: \n" << B << std::endl;
    std::cout   <<  "C: \n" << C << std::endl;
    std::cout   <<  "G: \n" << G << std::endl;

    //Observer obs(dt, A, B, C, G);
    Observer obs(dt, A, B, C, G);
    obs.init();
    

    // The following are dummy values

    Eigen::VectorXd u(30);
    Eigen::VectorXd Xc(30);
    Eigen::VectorXd Xz(30);

    std::cout << "u " << std::endl;

    u << 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 
         0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 
         0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1;

     std::cout << "Xc " << std::endl;

    Xc << 0.09, 0.214, 0.32, 0.382, 0.512, 0.584, 0.72, 0.818, 0.89, 1.02, 
          1.084, 1.19, 1.308, 1.404, 1.486, 1.592, 1.69, 1.808, 1.88, 2.004, 
          2.118, 2.192, 2.288, 2.408, 2.498, 2.6, 2.716, 2.788, 2.884, 3.012;

     std::cout << "Xz " << std::endl;

    Xz << 0.098, 0.192, 0.28, 0.386, 0.51, 0.582, 0.708, 0.818, 0.892, 0.994, 
          1.092, 1.196, 1.31, 1.398, 1.52, 1.586, 1.702, 1.786, 1.898, 2.01, 
          2.104, 2.2, 2.306, 2.398, 2.516, 2.58, 2.688, 2.802, 2.91, 2.982;
     

     for (size_t i = 0; i < 30; i++)
     {
          std::cout << "Timestamp: "  << i << std::endl;
          std::cout << "input u: "  << u[i] << std::endl;
          std::cout << "measure Xc: "  << Xc[i] << std::endl;
          std::cout << "measure Xz: "  << Xz[i] << std::endl;
          std::cout << "===================================" << std::endl;
     }

     std::cout << std::endl << std::endl;
     std::cout << "Init state: " << obs.state() << std::endl << std::endl;

     for (size_t i = 0; i < 30; i++)
     {
          Eigen::VectorXd u_i(1);
          Eigen::VectorXd y(2);
          u_i << u[i];
          y << Xc[i], Xz[i];

          std::cout << "Timestamp: "  << i << std::endl;
          std::cout << "Input u : "  << u_i << std::endl;
          std::cout << "Measures : "  << y << std::endl;
          obs.update(u_i, y);
          std::cout << "Updated State: " << obs.state() << std::endl;
          std::cout << "=====================================" << std::endl;
     }
     

}
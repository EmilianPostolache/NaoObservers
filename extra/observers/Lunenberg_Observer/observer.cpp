#include <iostream>
#include <stdexcept>

#include "observer.hpp"

Observer::Observer(double dt,
            const Eigen::MatrixXd& A,
            const Eigen::MatrixXd& B,
            const Eigen::MatrixXd& C,
            const Eigen::MatrixXd& G
        ): A(A), C(C), B(B), G(G), 
           m(C.rows()), n(A.rows()),initialized(false)
           {
               x_act.resize(n);
               x_next.resize(n);
           }

Observer::Observer() {}

void Observer::init(const Eigen::VectorXd& x0){
    x_act = x0;
    initialized = true;
}

void Observer::init(){
    x_act.setZero();
    initialized = true;
}

void Observer::update(const Eigen::VectorXd& u, const Eigen::VectorXd& y){
    if(!initialized)
        throw std::runtime_error("Filter is not initialized!");
    
    x_act = x_next;

    x_next = (A * x_act) + (B * u) + (G*(y - C*x_act));
}


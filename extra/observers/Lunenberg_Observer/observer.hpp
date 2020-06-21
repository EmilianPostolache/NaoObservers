#pragma once
#include "Eigen/Dense"

class Observer{
    public:
        Observer(double dt,
            const Eigen::MatrixXd& A,
            const Eigen::MatrixXd& B,
            const Eigen::MatrixXd& C,
            const Eigen::MatrixXd& G
        );

        Observer();

        void init();
        void init(const Eigen::VectorXd& x0);

        void update(const Eigen::VectorXd& u, const Eigen::VectorXd& y);

        Eigen::VectorXd state(){ return x_act;}


    private:
        // Matrices
        Eigen::MatrixXd A,B,C,G;

        // Estimated states 
        Eigen::VectorXd x_act, x_next;

        bool initialized;
        int m, n;
};
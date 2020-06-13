#pragma once
#include "Eigen/Dense"
#include <list>
#include <map>
#include <string>


class Observer{
    public:
        Observer(std::string name) { this->name = name; }
        virtual ~Observer() {}
        virtual void update(const Eigen::MatrixXd& U, const Eigen::MatrixXd& Y) = 0;
        virtual std::map<std::string, Eigen::VectorXd> state() = 0;
        std::string getName() { return name; }
    protected:
        std::string name;
};

class CompositeObserver : public Observer {
    public:
        CompositeObserver(std::string name): Observer(name) {};
        ~CompositeObserver();

        void add(Observer *obs);
        void rem(Observer *obs);

        void update(const Eigen::MatrixXd& U, const Eigen::MatrixXd& Y);
        std::map<std::string, Eigen::VectorXd> state();
    
    private:
        std::list<Observer *> children;
};

class LeafObserver : public Observer {
    public:
        LeafObserver(std::string name, int axis) :
                     Observer(name) {this->axis = axis;};
        virtual void init() { initialized = true; }

    protected:
        // Estimated states 
        int axis;
        bool initialized;
        Eigen::VectorXd xAct;
};


class LuenbergerObserver : public LeafObserver{
    public:
        LuenbergerObserver(
            const Eigen::MatrixXd& A,
            const Eigen::MatrixXd& B,
            const Eigen::MatrixXd& C,
            const Eigen::MatrixXd& G, std::string name, int axis);

        void init();
        void init(const Eigen::VectorXd& x0);
        void update(const Eigen::MatrixXd& U, const Eigen::MatrixXd& Y);
        std::map<std::string, Eigen::VectorXd> state();

     private:
         // Matrices
         Eigen::MatrixXd A,B,C,G;
         int m, n;
};

class KalmanFilter : public LeafObserver {

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

        KalmanFilter(const Eigen::MatrixXd& A,
                    const Eigen::MatrixXd& B,
                    const Eigen::MatrixXd& C,
                    const Eigen::MatrixXd& R,
                    double sigma_jerk,
                    double sigma_ddfext, std::string name, int axis);


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
        
        void update(const Eigen::MatrixXd& U, const Eigen::MatrixXd& Y);

        /*
        Update step (contains both prediction and correction) - for X, Y observers

        @param y: vector of measurements
        @param C: dynamic matrix to overwrite
        */

        void update(const Eigen::MatrixXd& C, const Eigen::MatrixXd& U, const Eigen::MatrixXd& Y);
        

        /*
        Returns the current estimate of the system state
        */
        std::map<std::string, Eigen::VectorXd> state();

        int getm(){ return m; }
        int getn(){ return n; }
        
    private:
        int m, n; // System dimensions
        Eigen::MatrixXd A, C, Q, R, P, K, NI; //Kalman Filter matrices
};

class KalmanComposite : public Observer {
    public:
        KalmanComposite(double g, double Mc): Observer("kalman_composite") {
            this->g = g;
            this->Mc = Mc;
        };
        ~KalmanComposite();

        void update(const Eigen::MatrixXd& U, const Eigen::MatrixXd& Y);
        std::map<std::string, Eigen::VectorXd> state();

        void addKalmanX(KalmanFilter *obs);
        void addKalmanY(KalmanFilter *obs);
        void addKalmanZ(KalmanFilter *obs);


    private:
        double g;
        double Mc;
        KalmanFilter *kalmanX;
        KalmanFilter *kalmanY;
        KalmanFilter *kalmanZ;
};


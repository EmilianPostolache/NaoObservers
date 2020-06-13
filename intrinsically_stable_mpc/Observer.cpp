#include <iostream>
#include <stdexcept>
#include "Observer.hpp"
#include <map>
#include <string>


/* ------------ Composite ---------------*/

void CompositeObserver::add(Observer *obs) {
     children.push_back(obs); }

CompositeObserver::~CompositeObserver () {
    for (Observer *c : children) {
        delete c;
    }
}

void CompositeObserver::rem(Observer *obs) { children.remove(obs); }

void CompositeObserver::update(const Eigen::MatrixXd& U, const Eigen::MatrixXd& Y) { 
    for (Observer *c : children) {
       c->update(U, Y);
    }
}

std::map<std::string, Eigen::VectorXd> CompositeObserver::state() {
    std::map<std::string, Eigen::VectorXd> m;
    for (Observer *c : children) {
        auto map = c->state();
        m.insert(map.begin(), map.end());
    }
    return m;
}

/* ---------- Lunenberg --------------- */

LuenbergerObserver::LuenbergerObserver(
            const Eigen::MatrixXd& A,
            const Eigen::MatrixXd& B,
            const Eigen::MatrixXd& C,
            const Eigen::MatrixXd& G, std::string name, int axis):
            LeafObserver(name, axis),
             A(A), C(C), B(B), G(G), 
           m(C.rows()), n(A.rows())
           {
               xAct.resize(n);
           }


void LuenbergerObserver::init(const Eigen::VectorXd& x0){
    LeafObserver::init();
    xAct = x0;
}

void LuenbergerObserver::init(){
    LeafObserver::init();
    xAct.setZero();
}

void LuenbergerObserver::update(const Eigen::MatrixXd& U, const Eigen::MatrixXd& Y){
    if(!initialized)
        throw std::runtime_error("Observer is not initialized!");
    Eigen::VectorXd u = U.row(this->axis);
    Eigen::Vector2d y;
    y << Y.row(this->axis)[0], Y.row(this->axis)[1];
    xAct = (A * xAct) + (B * u) + (G*(y - C*xAct));
}

std::map<std::string, Eigen::VectorXd> LuenbergerObserver::state(){
    std::map<std::string, Eigen::VectorXd> map;
    map[name] = xAct;
    return map;
}


/* ------------------ Kalman Composite ---------------- */


KalmanComposite::~KalmanComposite(){
    delete kalmanX;
    delete kalmanY;
    delete kalmanZ;
}


void KalmanComposite::update(const Eigen::MatrixXd& U, const Eigen::MatrixXd& Y) {
    kalmanZ->update(U, Y);

    double zcHat = kalmanZ->state()["kalman_z"](0);   // estimate of z CoM
    double ddzcHat = kalmanZ->state()["kalman_z"](2); // estimate of acceleration z CoM
    double fzHat = kalmanZ->state()["kalman_z"](3);   // estimate of external force z axis
    double grfHat = -Mc*g -Mc*ddzcHat + fzHat;

    Eigen::MatrixXd Cx(kalmanZ->getm(), kalmanZ->getn());
    Cx << 1, 0, 0, 0, 0,
          0, 0, 1, 0, 0,
          1, 0, Mc*zcHat/grfHat, -zcHat/grfHat, 0;

    kalmanX->update(Cx, U, Y);
    kalmanY->update(Cx, U, Y);
}

std::map<std::string, Eigen::VectorXd> KalmanComposite::state() {
    std::map<std::string, Eigen::VectorXd> m;
    auto mapX = kalmanX->state();
    m.insert(mapX.begin(), mapX.end());
    auto mapY = kalmanY->state();
    m.insert(mapY.begin(), mapY.end());
    auto mapZ = kalmanZ->state();
    m.insert(mapZ.begin(), mapZ.end());
    return m;
}
         
void KalmanComposite::addKalmanX(KalmanFilter *obs) {kalmanX = obs;}
void KalmanComposite::addKalmanY(KalmanFilter *obs) {kalmanY = obs;}
void KalmanComposite::addKalmanZ(KalmanFilter *obs) {kalmanZ = obs;}


/* ------------------ Kalman ------------------ */

// Constructor implementation
KalmanFilter::KalmanFilter(const Eigen::MatrixXd& A,
                            const Eigen::MatrixXd& B,
                            const Eigen::MatrixXd& C,
                            const Eigen::MatrixXd& R,
                            double sigma_jerk,
                            double sigma_ddfext, 
                            std::string name, int axis):LeafObserver(name, axis){

    Eigen::Matrix2d cov_input;
    cov_input << sigma_jerk, 0, sigma_ddfext, 0;
    Q = B*cov_input*B.transpose();                  // compute covariance process noise
    this->A = A;
    this->C = C;
    this->R = R;
    m = C.rows();
    n = A.rows();
    xAct.resize(n);
}


// KF initialization no guessing
void KalmanFilter::init(){

    int f = pow(10, 2);
    xAct.setZero();
    Eigen::MatrixXd p(n,n);
    P = p;
    P = f*P.setIdentity();
}

// KF initialization with guessing
void KalmanFilter::init(const Eigen::VectorXd& x0, const Eigen::MatrixXd& P0){
    xAct = x0;
    P = P0;
}

// Update function implementation (for Z observer)
void KalmanFilter::update(const Eigen::MatrixXd& U, const Eigen::MatrixXd& Y){
    Eigen::Vector3d y;
    y << Y.row(this->axis)[0], Y.row(this->axis)[2], Y.row(this->axis)[1];

    // Prediction step
    xAct = A*xAct;
    P = A*P*A.transpose() + Q;

    // Correction step
    K = P*C.transpose()*(C*P*C.transpose()+R).inverse();

    NI = y - C*xAct;
    xAct = xAct + K*NI;
    P = P - K*C*P;
}

// @overloaded update function implementation (for X, Y observers)
void KalmanFilter::update(const Eigen::MatrixXd& C, const Eigen::MatrixXd& U, const Eigen::MatrixXd& Y){
    this->C = C;
    update(U, Y);
}

// Returns current state estimate
std::map<std::string, Eigen::VectorXd> KalmanFilter::state(){
    std::map<std::string, Eigen::VectorXd> map;
    map[name] = xAct;
    return map;
}
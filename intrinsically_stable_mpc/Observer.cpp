#include <iostream>
#include <stdexcept>
#include "Observer.hpp"
#include <map>
#include <string>


void CompositeObserver::add(Observer *obs) { children.push_back(obs); }

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
               xNext.resize(n);
           }

LuenbergerObserver::LuenbergerObserver(std::string name, int axis)
             : LeafObserver(name, axis) {}


void LuenbergerObserver::init(const Eigen::VectorXd& x0){
    xAct = x0;
    LeafObserver::init(x0);
}

void LuenbergerObserver::init(){
    xAct.setZero();
    LeafObserver::init();
}

void LuenbergerObserver::update(const Eigen::MatrixXd& U, const Eigen::MatrixXd& Y){
    if(!initialized)
        throw std::runtime_error("Observer is not initialized!");
    Eigen::VectorXd u = U.row(this->axis);
    Eigen::VectorXd y = Y.row(this->axis);
    xNext = (A * xAct) + (B * u) + (G*(y - C*xAct));
    xAct = xNext;
}

std::map<std::string, Eigen::VectorXd> LuenbergerObserver::state(){
    std::map<std::string, Eigen::VectorXd> m;
    m[name] = xAct;
    return m;
}

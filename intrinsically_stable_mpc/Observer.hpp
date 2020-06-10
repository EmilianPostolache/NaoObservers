#pragma once
#include "Eigen/Dense"
#include <list>
#include <map>
#include <string>


class Observer{
    public:
        Observer(std::string name) { this->name = name; }
        Observer() {};
        virtual ~Observer() {}
        virtual void update(const Eigen::MatrixXd& u, const Eigen::MatrixXd& y) = 0;
        virtual std::map<std::string, Eigen::VectorXd> state() = 0;
    protected:
        std::string name;
};

class LeafObserver : public Observer {
    public:
        LeafObserver(std::string name, int axis) :
                     Observer(name) {this->axis = axis;};
        LeafObserver() {};  
        virtual void init() { initialized = true; }
        virtual void init(const Eigen::VectorXd& x0){ initialized = true; }

    protected:
        // Estimated states 
        int axis;
        bool initialized;
        Eigen::VectorXd xAct, xNext;
};

class CompositeObserver : public Observer {
    public:
        CompositeObserver(std::string name): Observer(name) {};
        CompositeObserver() {};

        void add(Observer *obs);
        void rem(Observer *obs);

        void update(const Eigen::MatrixXd& u, const Eigen::MatrixXd& y);
        std::map<std::string, Eigen::VectorXd> state();
    
    private:
        std::list<Observer *> children;
};

class LuenbergerObserver : public LeafObserver{
    public:
        LuenbergerObserver(
            const Eigen::MatrixXd& A,
            const Eigen::MatrixXd& B,
            const Eigen::MatrixXd& C,
            const Eigen::MatrixXd& G, std::string name, int axis);
        LuenbergerObserver(std::string name, int axis);
        LuenbergerObserver() {};
        void init();
        void init(const Eigen::VectorXd& x0);
        void update(const Eigen::MatrixXd& u, const Eigen::MatrixXd& y);
        std::map<std::string, Eigen::VectorXd> state();

     private:
         // Matrices
         Eigen::MatrixXd A,B,C,G;
         int m, n;
};

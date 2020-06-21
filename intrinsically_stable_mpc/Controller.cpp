#include "Controller.hpp"
#include "utils.cpp"
#include "Observer.hpp"
#include <iostream>
#include <stdlib.h>
#include <cmath>
#include <map>

using namespace std;

double body_angles[] = {0,0,0};

//TODO: loggare le covarianze (P)
//TODO: sistemare il vettore delle forze nei log
//TODO: implementare il terzo observer

Controller::Controller(dart::dynamics::SkeletonPtr _robot,
		dart::simulation::WorldPtr _world)
: mRobot(_robot),
  mWorld(_world)
{
	assert(_robot != nullptr);

	ikGain = 20;

	// some useful pointers to robot limbs
	leftFoot = mRobot->getBodyNode("l_sole");
	rightFoot = mRobot->getBodyNode("r_sole");
	mBase = mRobot->getBodyNode("base_link");
	mTorso = mRobot->getBodyNode("torso");

	supportFoot = LEFT;
	swingFootStartingPosition.resize(6);
	if (supportFoot == RIGHT) {
		mSupportFoot = mRobot->getBodyNode("r_sole");
		mSwingFoot = mRobot->getBodyNode("l_sole");
	} else {
		mSupportFoot = mRobot->getBodyNode("l_sole");
		mSwingFoot = mRobot->getBodyNode("r_sole");
	}

	balancePoint = COM;
    beheavior = WALK;

	singleSupportDuration = 0.3;
	doubleSupportDuration = 0.2;

	stepHeight = 0.02;

	indInitial = mWorld->getSimFrames() + (int)((singleSupportDuration + doubleSupportDuration)/mWorld->getTimeStep());

	Eigen::Vector3d comInitialPosition;
	if (balancePoint == TORSO) {
		comInitialPosition = mBase->getCOM(mSupportFoot);
	} else {
		comInitialPosition = mRobot->getCOM(mSupportFoot);
	}
	lastMeasZmp << mTorso->getCOM()[0], mTorso->getCOM()[1], 0.0;
	measZmp = lastMeasZmp;

	double comTargetHeight = 0.33;//0.33;

	solver = new mpcSolver::MPCSolver(0.05, mWorld->getTimeStep(), 1.0, comInitialPosition, comTargetHeight, singleSupportDuration, doubleSupportDuration, 0.30,
			0.05, 0.2, 0.1, 0.15, 0.1, 0.0); //0.05, 0.2, 0.1, 0.15, 0.1, 0.0   0.05, 0.3, 0.13, 0.17, 0.1, 0.0

	swingFootStartingPosition << getRPY(mSwingFoot, mSupportFoot), mSwingFoot->getCOM(mSupportFoot);

	//Balance
	balanceBasePos.resize(6);
	balanceBasePos << getRPY(mBase, mSupportFoot), comInitialPosition;
	balanceBasePos(5) = comTargetHeight;
	balanceFootPos.resize(6);
	balanceFootPos = swingFootStartingPosition;

	setObservers(comTargetHeight, comInitialPosition, 1.0/100.0);
}

Controller::~Controller()
{
	delete observers;
}

void Controller::setObservers(double comTargetHeight, Eigen::Vector3d& comInitialPosition, double dt){
	
	// ------------------------- Observers --------------------------

	observers = new CompositeObserver("observers");
    double ni = sqrt(g/comTargetHeight);

	//------------------------ Luenberger Observer ------------------------

	int nLun = 5;
    int mLun = 2;

    Eigen::MatrixXd ALun(nLun, nLun);
    Eigen::MatrixXd BLun(nLun, 1);
    Eigen::MatrixXd CLun(mLun, nLun);
    Eigen::MatrixXd GLun(nLun, mLun);


    // ALun << 1, dt, 0, 0, 0, 
    //         pow(ni,2)*dt, 1, -pow(ni,2)*dt, dt, 0,
    //         0, 0, 1, 0, 0,
    //         0, 0, 0, 1, dt,
    //         0, 0, 0, 0, 1;

    // BLun << 0, 0, dt, 0, 0;

    // CLun << 1, 0, 0, 0, 0,
    //         0, 0, 1, 0, 0;

    // GLun << -2.2528, -2.5477,
	// 		31.4317, -25.7226,
	// 	    -0.0001, -0.7972,
	// 	    -0.4593, -1.6285,
	// 		 0.0190, 0.0722;

	ALun << 0, 1, 0, 0, 0,
	        pow(ni,2), 0, -pow(ni,2), 1, 0,
			0, 0, 0, 0, 0,
			0, 0, 0, 0, 1,
			0, 0, 0, 0, 0;

	BLun << 0, 0, 1, 0, 0;

	CLun << 1, 0, 0, 0, 0,
	        0, 0, 1, 0, 0;


    // GLun << 158.82, -444.63,
	// 		9371.79, -51423.15,
	// 		-0.02, 41.18,
	// 		241151.90, -1941421.30,
	// 		2304743.72, -24002899.99;

	GLun << 330.00, 0.00,
			40804.70, -29.70,
			-0.00, 70.00,
			2235750.00, 1.00,
			45900000.14, 70.00;


    LuenbergerObserver *xLunObs = new LuenbergerObserver(ALun, BLun, CLun, GLun, "luenberger_x", 0);
    LuenbergerObserver *yLunObs = new LuenbergerObserver(ALun, BLun, CLun, GLun, "luenberger_y", 1);
	LuenbergerObserver *zLunObs = new LuenbergerObserver(ALun, BLun, CLun, GLun, "luenberger_z", 2);

    observers->add(xLunObs);
    observers->add(yLunObs);
	observers->add(zLunObs);
	// observers.add(lzObs);

    Eigen::VectorXd initXLun (nLun);
	Eigen::VectorXd initYLun (nLun);
	Eigen::VectorXd initZLun (nLun);

	// Eigen::VectorXd initZ (n);
    initXLun << comInitialPosition[0], 0., comInitialPosition[0], 0., 0.;
    initYLun << comInitialPosition[1], 0., comInitialPosition[1], 0., 0.;
	initZLun << comInitialPosition[2], 0., 0., 0., 0.;
	
	xLunObs -> init(initXLun);
 	yLunObs -> init(initYLun);
	zLunObs -> init(initZLun);

	//------------------------ Kalman Observer ------------------------

	int uKal = 2;          // size input vector
	int nKal = 5;          // size vector states
	int mKal = 3;          // size vector measurements

	double sigmaQJerk = pow(10.0, 3.0);
	double sigmaQDdfext = pow(10.0, 3.0);

	double zcHat, ddzcHat, fzHat, grfHat;
	double fxHat, fyHat, ts;

	// Initialize matrices for observers

	Eigen::MatrixXd RzKal(mKal, mKal);
	Eigen::MatrixXd RxKal(mKal, mKal);
	Eigen::MatrixXd RyKal(mKal, mKal);

	Eigen::MatrixXd AKal(nKal, nKal);
	Eigen::MatrixXd BKal(nKal, uKal);
	Eigen::MatrixXd CzKal(mKal, nKal);

	// Define matrices for observers

	RzKal << 0.01, 0, 0, 
	          0, 1, 0,
			  0, 0, 1;
	RxKal << 0.01, 0, 0,
	          0, 1, 0,
			  0, 0, 0.01;
	RyKal = RxKal;

	AKal << 1, dt, pow(dt, 2.0)/2, 0, 0,
            0, 1, dt, 0, 0, 
            0, 0, 1, 0, 0, 
            0, 0, 0, 1, dt, 
            0, 0, 0, 0, 1;
            
    BKal << pow(dt, 3.0)/6, 0,
		   pow(dt, 2.0)/2, 0,
		   dt, 0, 
		   0, pow(dt, 2.0)/2,
		   0, dt;

    CzKal << 1, 0, 0, 0, 0, 
            0, 0, 1, 0, 0, 
            0, 0, -Mc, 1, 0;

	// dargli delle inizializzazioni
	KalmanFilter* zKalObs = new KalmanFilter(AKal, BKal, CzKal, RzKal, sigmaQJerk,
	                      sigmaQDdfext, "kalman_z", 2);  // z-axis Observer
	zKalObs->init();

	KalmanFilter* xKalObs = new KalmanFilter(AKal, BKal, CzKal, RxKal, sigmaQJerk,
	                     sigmaQDdfext, "kalman_x", 0);  // x-axis Observer
	xKalObs->init();

	KalmanFilter* yKalObs = new KalmanFilter(AKal, BKal, CzKal, RyKal, sigmaQJerk,
	                     sigmaQDdfext, "kalman_y", 1);  // z-axis Observer
	yKalObs->init();
	KalmanComposite* kalObs = new KalmanComposite(g, Mc);
	kalObs->addKalmanX(xKalObs);
	kalObs->addKalmanY(yKalObs);
	kalObs->addKalmanZ(zKalObs);

	observers->add(kalObs);

	/* --------------------- Stephens Filter ------------------ */

	int nSt = 4;          // size vector states
    int mSt = 2;          // size vector measurements
    int uSt = 1;          // size input vector

    double posProNoise = exp(-8), 
	velProNoise = exp(-4);
	double posOutputNoise = exp(-5),
	zmpOutputNoise = exp(-4);

    // Initialize matrices for observers
    Eigen::MatrixXd Q(nSt, nSt);
    Eigen::MatrixXd R(mSt, mSt);

    Eigen::MatrixXd A(nSt, nSt);
    Eigen::MatrixXd B(nSt, uSt);
    Eigen::MatrixXd C(mSt, nSt);

    // Define matrices for observers
    Q <<  posProNoise, 0, 0, 0,
		  0, velProNoise, 0, 0,
          0, 0, velProNoise, 0,
          0, 0, 0, velProNoise;

    R << posOutputNoise, 0,
		 0, posOutputNoise;

    A << 1,           dt,             0, 0, 
         pow(ni,2)*dt, 1, -pow(ni,2)*dt, dt,
		 0,            0,             1, 0, 
         0,            0,             0, 1;

        
    B << 0, 0, dt, 0;

    C << 1, 0, 0, 0, 
		 0, 0, 1, 0;


	int f = pow(10, 2);
	Eigen::MatrixXd P(nSt, nSt);
    P = f*P.setIdentity();
	Eigen::VectorXd initXSt (nSt);
	Eigen::VectorXd initYSt (nSt);
	Eigen::VectorXd initZSt (nSt);

	// Eigen::VectorXd initZ (n);
    initXSt << comInitialPosition[0], 0., comInitialPosition[0], 0.;
    initYSt << comInitialPosition[1], 0., comInitialPosition[1], 0.;
	initZSt << comInitialPosition[2], 0., 0., 0.;

	// ---------- rest 

    // Initialization Observers and vector of measurements
    StephensFilter *xStObs = new StephensFilter(A, B, C, Q, R, "stephens_x", 0); 
    xStObs->init(initXSt, P);
	StephensFilter *yStObs = new StephensFilter(A, B, C, Q, R, "stephens_y", 1);  
    yStObs->init(initYSt, P);
  	StephensFilter *zStObs = new StephensFilter(A, B, C, Q, R, "stephens_z", 2);  
    zStObs->init(initZSt, P);

	observers->add(xStObs);
	observers->add(yStObs);
	observers->add(zStObs);
}

void Controller::setExternalForceConstant(){
	externalForceMode = 0;
}

void Controller::setExternalForcePeriodic(){
	externalForceMode = 1;
}

void Controller::setExternalForceStartFrame(int startFrame){
	externalForceStartFrame = startFrame;
}

void Controller::setExternalForceX(float f){
	externalForce[0] = f;
}

void Controller::setExternalForceY(float f){
	externalForce[1] = f;
}

void Controller::setExternalForceZ(float f){
	externalForce[2] = f;
}

void Controller::setExternalForcePeriodicPhase(float phi){
	externalForcePeriodicPhase = phi;
}

void Controller::setExternalForcePeriodicFrequency(float k){
	externalForcePeriodicFrequency = k;
}

Eigen::Vector3d Controller::getExternalForce(){
	if (mWorld->getSimFrames()>=externalForceStartFrame && mWorld->getSimFrames()<=10000000) {
		switch (externalForceMode){
			case 0:
				// Constant force
				return externalForce;
			case 1:
				// Periodic force
				return Eigen::Vector3d(externalForce[0]*cos(2*M_PI * externalForcePeriodicFrequency * mWorld->getSimFrames()*0.01 + externalForcePeriodicPhase),
									externalForce[1]*sin(2*M_PI * externalForcePeriodicFrequency * mWorld->getSimFrames()*0.01 + externalForcePeriodicPhase), 0.0);
			default:
				return Eigen::Vector3d(0.0, 0.0, 0.0);
      }
	}
	return Eigen::Vector3d(0.0, 0.0, 0.0);
}

void Controller::update()
{ 
    std::random_device dev;
    std::mt19937 rng(dev());
    std::uniform_int_distribution<std::mt19937::result_type> dist6(0,15); // distribution in range [1, 6]

    // std::cout << "AAA"<<mWorld->getSimFrames() << std::endl;
	Eigen::MatrixXd Y (3, 3); 
	Eigen::MatrixXd U (3, 1);
	// std::pair<Eigen::Vector3d, double> wrench = this->getZmpFromWrench();
	// double grf = wrench.second;
	//Eigen::Vector3d zmpMeas = wrench.first;
	std::pair<Eigen::Vector3d, Eigen::Vector2d> fromExtForces = this->getZmpFromExternalForces();
	measZmp = fromExtForces.first;
	Eigen::Vector2d grf =fromExtForces.second;
	// TODO: accertarsi che non andra' in negativo
	if (measZmp[0] == 0) {
		measZmp[0] = lastMeasZmp[0];
	}
	if (measZmp[1] == 0) {
		measZmp[1] = lastMeasZmp[1];
	}
	lastMeasZmp = measZmp;


    //Eigen::Vector3d zmpMeas = this->getZmpFromAngularMomentum();

	//Eigen::Vector3d zmpMeas = this->getZmpFromExternalForces();


	// zmpMeas = mSupportFoot->getCOM() + (mSupportFoot->getWorldTransform().rotation())*(solver->getOptimalZMPPosition());
	
	U = solver -> getZMPDot();
	Eigen::Vector3d comCurrentPosition;
	Eigen::VectorXd comCurrentAcceleration;
	Eigen::VectorXd comCurrentVelocity;

	// if (balancePoint == TORSO) {
	//   comCurrentPosition = mBase->getCOM(/*mSupportFoot*/) - mSupportFoot->getCOM();
	//	 comCurrentAcceleration = mBase->getCOMLinearAcceleration(/*mSupportFoot*/);
	// } 
	// else {
	//    comCurrentPosition = mRobot->getCOM(/*mSupportFoot*/) - mSupportFoot->getCOM();
	//	comCurrentAcceleration = mRobot->getCOMLinearAcceleration(/*mSupportFoot*/);
	// }
 	comCurrentPosition = mTorso->getCOM(); //- mSupportFoot->getCOM();
	comCurrentAcceleration = mTorso->getCOMLinearAcceleration();
	comCurrentVelocity = mTorso->getCOMLinearVelocity();
	// zmpMeas -= mSupportFoot->getCOM();
	Y << comCurrentPosition[0], measZmp[0], comCurrentAcceleration[0],
         comCurrentPosition[1], measZmp[1], comCurrentAcceleration[1],
	 	 comCurrentPosition[2], -(grf[0]+grf[1]), comCurrentAcceleration[2];

	observers->update(U, Y);
	mTorso->addExtForce(getExternalForce());//, mTorso->getLocalCOM());

	// if (mWorld->getSimFrames()>=250 && mWorld->getSimFrames()<=10000000) 
	// 	mTorso->addExtForce(Eigen::Vector3d(290*cos(2*3.14*((mWorld->getSimFrames())*6.28)/3),290*sin(2*3.14*((mWorld->getSimFrames())*3.14)/0.3),0));  //0.86*sin(2*3.14*(mWorld->getTimeStep())/0.2)*180 magic value!
    // //if (mWorld->getSimFrames()>=250 && mWorld->getSimFrames()<=10000000) mSupportFoot->addExtForce(Eigen::Vector3d(10-15*cos(2*3.14*((mWorld->getSimFrames())*6.28)/30),-10-15*sin(2*3.14*((mWorld->getSimFrames())*3.14)/40),0));  //150 70
	// if (mWorld->getSimFrames()>=250 && mWorld->getSimFrames()<=10000000) mTorso->addExtForce(Eigen::Vector3d(-15+10*cos(2*3.14*((mWorld->getSimFrames())*6.28)/3),35+15*sin(2*3.14*((mWorld->getSimFrames())*3.14)/0.3),0)); 
    
    // # sync; echo 3 > /proc/sys/vm/drop_cache
	Eigen::VectorXd qDot;
	if (beheavior == WALK) qDot = generateWalking();
	else qDot = generateBalance();
	
	// TODO: where is the point in which data has to be used?? Two possible cases,
	// here or after setCommand has been issued


	// Set the velocity of the floating base to zero
	for (int i=0; i<6; ++i){
		mRobot->setCommand(i, 0);
	}

	// Set the velocity of each joint as per inverse kinematics
	for (int i=0; i<24; ++i){
		mRobot->setCommand(i+6,qDot(i));
	}

	qDotOld = qDot;

	// Store the results in files (for plotting)
	storeData();
}

void Controller::storeData() {

	Eigen::Vector3d worldPredictedCom = mSupportFoot->getCOM() + (mSupportFoot->getWorldTransform().rotation())*(solver->getOptimalCoMPosition());
	Eigen::Vector3d worldPredictedZmp = mSupportFoot->getCOM() + (mSupportFoot->getWorldTransform().rotation())*(solver->getOptimalZMPPosition());

	double filterGain = 0.2;
	//zmpPosMeas = mSupportFoot->getCOM() + (mSupportFoot->getWorldTransform().rotation())*(zmpPosMeas);

	// std::ofstream fout0("/home/emilian/Desktop/Mobile Robotics/project/NaoObservers/intrinsically_stable_mpc/data/zmp_meas_x.txt", std::ofstream::app);
	// fout0 << measZmp[0] << std::endl;

	// std::ofstream fout1("/home/emilian/Desktop/Mobile Robotics/project/NaoObservers/intrinsically_stable_mpc/data/zmp_meas_y.txt", std::ofstream::app);
	// fout1 << measZmp[1] << std::endl;

	//std::ofstream fout3("/home/emilian/Desktop/Mobile Robotics/project/NaoObservers/intrinsically_stable_mpc/data/comPositionNominal.txt", std::ofstream::app);
   	//fout3 << (mSupportFoot->getCOM()+solver->getOptimalCoMPosition()).transpose()/*worldPredictedCom.transpose()*/ << std::endl;

	// std::ofstream fout4("/home/emilian/Desktop/Mobile Robotics/project/NaoObservers/intrinsically_stable_mpc/data/comPositionMeasured.txt", std::ofstream::app);
   	// fout4 << (mSupportFoot->getCOM() + mRobot->getCOM()).transpose() << std::endl;

	// std::ofstream fout5("/home/emilian/Desktop/Mobile Robotics/project/NaoObservers/intrinsically_stable_mpc/data/comVelocityNominal.txt", std::ofstream::app);
   	// fout5 << solver->getOptimalCoMVelocity().transpose() << std::endl;

	// std::ofstream fout6("/home/emilian/Desktop/Mobile Robotics/project/NaoObservers/intrinsically_stable_mpc/data/comVelocityMeasured.txt", std::ofstream::app);
   	// fout6 << mRobot->getCOMLinearVelocity().transpose() << std::endl;

	// std::ofstream fout7("/home/emilian/Desktop/Mobile Robotics/project/NaoObservers/intrinsically_stable_mpc/data/zmpPositionNominal.txt", std::ofstream::app);
   	// fout7 << (mSupportFoot->getCOM()+solver->getOptimalZMPPosition()).transpose()/*worldPredictedZmp.transpose()*/ << std::endl;

	// std::ofstream fout8("/home/emilian/Desktop/Mobile Robotics/project/NaoObservers/intrinsically_stable_mpc/data/zmpPositionMeasured.txt", std::ofstream::app);
   	// fout8 << (mSupportFoot->getCOM() + mRobot->getCOM() - mSupportFoot->getCOM()).transpose() - mRobot->getCOMLinearAcceleration().transpose()/(solver->getOmega()*solver->getOmega()) << std::endl;


	std::map<std::string, Eigen::VectorXd> states = observers->state();
	std::map<std::string, Eigen::MatrixXd> uncertainty = observers->uncertainty();

	//---------------------- luenberger ----------------------
	Eigen::Vector3d force = getExternalForce();


	std::ofstream fout9("/home/emilian/Desktop/Mobile Robotics/project/NaoObservers/intrinsically_stable_mpc/data/luenberger_x.txt", std::ofstream::app);
	Eigen::VectorXd state = states["luenberger_x"];
	state[3] *= Mc;
	fout9 << state.transpose() << std::endl;

	std::ofstream fout10("/home/emilian/Desktop/Mobile Robotics/project/NaoObservers/intrinsically_stable_mpc/data/luenberger_x_gt.txt", std::ofstream::app);
	fout10 << mTorso->getCOM()[0] << " " //- mSupportFoot->getCOM()[0] << " "
	       << mTorso->getCOMLinearVelocity()[0] << " " 
		   << measZmp[0] <<  " "
		   << force[0] <<  " "
		   << 0.0 <<std::endl;

	std::ofstream fout11("/home/emilian/Desktop/Mobile Robotics/project/NaoObservers/intrinsically_stable_mpc/data/luenberger_y.txt", std::ofstream::app);
	// fout11 << states["luenberger_y"].transpose() << std::endl;
	state = states["luenberger_y"];
	state[3] *= Mc;
	fout11 << state.transpose() << std::endl;

	std::ofstream fout12("/home/emilian/Desktop/Mobile Robotics/project/NaoObservers/intrinsically_stable_mpc/data/luenberger_y_gt.txt", std::ofstream::app);
	fout12 << mTorso->getCOM()[1] << " " // - mSupportFoot->getCOM()[1] 
	       << mTorso->getCOMLinearVelocity()[1] << " " 
		   << measZmp[1] <<  " "
		   << force[1] <<  " "
		   << 0.0 <<std::endl;


	// std::ofstream fout13("/home/emilian/Desktop/Mobile Robotics/project/NaoObservers/intrinsically_stable_mpc/data/luenberger_z.txt", std::ofstream::app);
	// fout13 << states["luenberger_z"].transpose() << std::endl;

	// std::ofstream fout14("/home/emilian/Desktop/Mobile Robotics/project/NaoObservers/intrinsically_stable_mpc/data/luenberger_z_gt.txt", std::ofstream::app);
	// fout14 << mTorso->getCOM()[2] << " " // - mSupportFoot->getCOM()[2] << " "
	//        << mTorso->getCOMLinearVelocity()[2] << " " 
	// 	   << measZmp[2] <<  " "
	// 	   << force[2] <<  " "
	// 	   << 0.0 <<std::endl;

	//---------------------- kalman ----------------------

	//TODO: non scordarci di aggiungere le derivate corrette della forza misurate

	std::ofstream fout15("/home/emilian/Desktop/Mobile Robotics/project/NaoObservers/intrinsically_stable_mpc/data/kalman_x.txt", std::ofstream::app);
	fout15 << states["kalman_x"].transpose() << std::endl;

	std::ofstream fout16("/home/emilian/Desktop/Mobile Robotics/project/NaoObservers/intrinsically_stable_mpc/data/kalman_x_gt.txt", std::ofstream::app);
	fout16 << mTorso->getCOM()[0] << " " //- mSupportFoot->getCOM()[0] << " "
	       << mTorso->getCOMLinearVelocity()[0] << " " 
		   << mTorso->getCOMLinearAcceleration()[0] << " "
		   << force[0] <<  " "
		   << measZmp[0] << std::endl;

	std::ofstream fout166("/home/emilian/Desktop/Mobile Robotics/project/NaoObservers/intrinsically_stable_mpc/data/kalman_x_un.txt", std::ofstream::app);
	fout166 << uncertainty["kalman_x"] << std::endl;


	std::ofstream fout17("/home/emilian/Desktop/Mobile Robotics/project/NaoObservers/intrinsically_stable_mpc/data/kalman_y.txt", std::ofstream::app);
	fout17 << states["kalman_y"].transpose() << std::endl;
	
	std::ofstream fout18("/home/emilian/Desktop/Mobile Robotics/project/NaoObservers/intrinsically_stable_mpc/data/kalman_y_gt.txt", std::ofstream::app);
	fout18 << mTorso->getCOM()[1] << " " // - mSupportFoot->getCOM()[1] << " "
		<< mTorso->getCOMLinearVelocity()[1] << " " 
		<< mTorso->getCOMLinearAcceleration()[1] << " "
		<< force[1] <<  " "
		<< measZmp[1] << std::endl;

	std::ofstream fout188("/home/emilian/Desktop/Mobile Robotics/project/NaoObservers/intrinsically_stable_mpc/data/kalman_y_un.txt", std::ofstream::app);
	fout188 << uncertainty["kalman_y"] << std::endl;

	std::ofstream fout19("/home/emilian/Desktop/Mobile Robotics/project/NaoObservers/intrinsically_stable_mpc/data/kalman_z.txt", std::ofstream::app);

	Eigen::VectorXd stateKalmanZ = states["kalman_z"];
	//stateKalmanZ[4] = -Mc*g - Mc * stateKalmanZ[2] + stateKalmanZ[3];
	fout19 << stateKalmanZ.transpose() << std::endl;

	//std::pair<Eigen::Vector3d, double> wrench = this->getZmpFromWrench();
	//double grf = -wrench.second;
	std::pair<Eigen::Vector3d, Eigen::Vector2d> fromExtForces = this->getZmpFromExternalForces();
	Eigen::Vector2d measGrf = fromExtForces.second;

	std::ofstream fout20("/home/emilian/Desktop/Mobile Robotics/project/NaoObservers/intrinsically_stable_mpc/data/kalman_z_gt.txt", std::ofstream::app);
	fout20 << mTorso->getCOM()[2] << " " // - mSupportFoot->getCOM()[2] << " "
		<< mTorso->getCOMLinearVelocity()[2] << " " 
		<< mTorso->getCOMLinearAcceleration()[2] << " "
		<< force[2] <<  " "
		<< -(measGrf[0]+measGrf[1])*Mc+Mc*g << std::endl;
	
	std::ofstream fout200("/home/emilian/Desktop/Mobile Robotics/project/NaoObservers/intrinsically_stable_mpc/data/kalman_z_un.txt", std::ofstream::app);
	fout200 << uncertainty["kalman_z"] << std::endl;

	//---------------------- stephens ----------------------

	std::ofstream fout21("/home/emilian/Desktop/Mobile Robotics/project/NaoObservers/intrinsically_stable_mpc/data/stephens_x.txt", std::ofstream::app);
	state = states["stephens_x"];
	state[3] *= Mc;
	fout21 << state.transpose() << std::endl;

	std::ofstream fout22("/home/emilian/Desktop/Mobile Robotics/project/NaoObservers/intrinsically_stable_mpc/data/stephens_x_gt.txt", std::ofstream::app);
	fout22 << mTorso->getCOM()[0] <<  " " //- mSupportFoot->getCOM()[0] << " "
	       << mTorso->getCOMLinearVelocity(/*mSupportFoot*/)[0] << " " 
		   << measZmp[0] <<  " " //- mSupportFoot->getCOM()[0] <<  " "
		   << force[0] << std::endl;

	std::ofstream fout222("/home/emilian/Desktop/Mobile Robotics/project/NaoObservers/intrinsically_stable_mpc/data/stephens_x_un.txt", std::ofstream::app);
	fout222 << uncertainty["stephens_x"] << std::endl;

	std::ofstream fout23("/home/emilian/Desktop/Mobile Robotics/project/NaoObservers/intrinsically_stable_mpc/data/stephens_y.txt", std::ofstream::app);
	state = states["stephens_y"];
	state[3] *= Mc;
	fout23 << state.transpose() << std::endl;
	
	std::ofstream fout24("/home/emilian/Desktop/Mobile Robotics/project/NaoObservers/intrinsically_stable_mpc/data/stephens_y_gt.txt", std::ofstream::app);
	fout24 << mTorso->getCOM()[1]  <<  " " //- mSupportFoot->getCOM()[1] << " "
	       << mTorso->getCOMLinearVelocity(/*mSupportFoot*/)[1] << " " 
		   << measZmp[1] <<  " " //- mSupportFoot->getCOM()[1]<<  " "
		   << force[1] << std::endl;

	std::ofstream fout244("/home/emilian/Desktop/Mobile Robotics/project/NaoObservers/intrinsically_stable_mpc/data/stephens_y_un.txt", std::ofstream::app);
	fout244 << uncertainty["stephens_y"] << std::endl;

	// std::ofstream fout25("/home/emilian/Desktop/Mobile Robotics/project/NaoObservers/intrinsically_stable_mpc/data/stephens_z.txt", std::ofstream::app);
	// fout25 << states["stephens_z"].transpose() << std::endl;

	// std::ofstream fout26("/home/emilian/Desktop/Mobile Robotics/project/NaoObservers/intrinsically_stable_mpc/data/stephens_z_gt.txt", std::ofstream::app);
	// fout26 << mTorso->getCOM()[2] <<  " " //- mSupportFoot->getCOM()[2] << " "
	//        << mTorso->getCOMLinearVelocity()[2] << " " 
	// 	   << measZmp[2]  <<  " " //- mSupportFoot->getCOM()[2] <<  " "
	// 	   << externalForce[2] << std::endl;
}

Eigen::VectorXd Controller::generateWalking(){

	if (solver->supportFootHasChanged()) supportFoot = !supportFoot;

	if (supportFoot == RIGHT) {
		mSupportFoot = mRobot->getBodyNode("r_sole");
		mSwingFoot = mRobot->getBodyNode("l_sole");
	} else {
		mSupportFoot = mRobot->getBodyNode("l_sole");
		mSwingFoot = mRobot->getBodyNode("r_sole");
	}

	if (solver->supportFootHasChanged()) {
		swingFootStartingPosition << getRPY(mSwingFoot, mSupportFoot), mSwingFoot->getCOM(mSupportFoot);
		if (mWorld->getSimFrames()>0) indInitial = mWorld->getSimFrames();
		footstepCounter++;
	}

	// Retrieve positions and orientations relative to support foot
	Eigen::VectorXd comCurrentPosition;
	Eigen::VectorXd comCurrentVelocity;
	Eigen::VectorXd comCurrentAcceleration;
	if (balancePoint == TORSO) {
		comCurrentPosition = mBase->getCOM(/*mSupportFoot*/) - mSupportFoot->getCOM();
		comCurrentVelocity = mBase->getCOMLinearVelocity(/*mSupportFoot*/);
		comCurrentAcceleration = mBase->getCOMLinearAcceleration(/*mSupportFoot*/);
	} else {
		comCurrentPosition = mRobot->getCOM(/*mSupportFoot*/) - mSupportFoot->getCOM();
		comCurrentVelocity = mRobot->getCOMLinearVelocity(/*mSupportFoot*/);
		comCurrentAcceleration = mRobot->getCOMLinearAcceleration(/*mSupportFoot*/);
	}

  std::random_device dev;
    std::mt19937 rng(dev());
    std::uniform_int_distribution<std::mt19937::result_type> dist6(-15,15); // distribution in range [1, 6]
   // +0.15*(dist6(rng)/15)
   //      	std::cout << "comCurrentP " << comCurrentPosition(2) <<"time "<< mWorld->getSimFrames()<< std::endl;

   //if (mWorld->getSimFrames()>=500 && mWorld->getSimFrames()<=50020) comCurrentPosition(0) = comCurrentPosition(0) + 0*0.0001*sin(mWorld->getSimFrames()-500);
   //if (mWorld->getSimFrames()>=500 && mWorld->getSimFrames()<=50020) comCurrentPosition(1) = comCurrentPosition(1) + 0*0.0001*cos((mWorld->getSimFrames())/0.3);
   //if (mWorld->getSimFrames()>=500 && mWorld->getSimFrames()<=50020) comCurrentPosition(2) = comCurrentPosition(2) - 0*0.0001*cos((mWorld->getSimFrames())/0.3);

   //     	std::cout << "comCurrentP plus " << comCurrentPosition(2) << std::endl;


	Eigen::VectorXd actPosSwingFoot(6);
	Eigen::VectorXd actPosBase(6);
	actPosSwingFoot << getRPY(mSwingFoot, mSupportFoot), mSwingFoot->getCOM(mSupportFoot);
	actPosBase << getRPY(mBase, mSupportFoot), comCurrentPosition;

	Eigen::Vector3d zmpCurrentPosition = Eigen::Vector3d::Zero();
	/*if (mWorld->getSimFrames() > 0) {
		zmpCurrentPosition = Eigen::Vector3d(zmpPosMeasStoreFiltered[0].back(), zmpPosMeasStoreFiltered[1].back(), 0);
		zmpCurrentPosition = -mSupportFoot->getCOM() + (mSupportFoot->getWorldTransform().rotation().transpose())*zmpCurrentPosition;
	}*/

	// Compute the CoM prediction using MPC
	solver->solve(comCurrentPosition, comCurrentVelocity, comCurrentAcceleration, mSwingFoot->getTransform(mSupportFoot),
			supportFoot, mWorld->getTime(), 0.1, 0.0, 0.0); //max vel is 0.16, max omega is 0.1
	// std::cout << "comCurrentPos" << comCurrentPosition << std::endl;
	// std::cout << "comCurrentVelocity" << comCurrentVelocity << std::endl;
	// std::cout << "comCurrentAcceleration" << comCurrentAcceleration << std::endl;
	Eigen::Affine3d temp = mSwingFoot->getTransform(mSupportFoot);
	// std::cout << "mSwing" << temp << std::endl;
	// std::cout << "supportFoot" << supportFoot << std::endl;
	// Assemble desired tasks
	Eigen::VectorXd desPosSwingFoot = getOmnidirectionalSwingFootTrajectoryMPC(solver->getOptimalFootsteps(),swingFootStartingPosition, actPosSwingFoot, stepHeight, indInitial);

        Eigen::VectorXd OptComPos(3);
        OptComPos << solver->getOptimalCoMPosition();


        // Data for angular momentum feedback

        Eigen::VectorXd vel(3);
        vel << mBase->getAngularVelocity();
        body_angles[0] = body_angles[0]+0.05*vel[0];
        body_angles[1] = body_angles[1]+0.05*vel[1];
        body_angles[2] = body_angles[2]+0.05*vel[2];

    //    std::cout << "Body Angles" << std::endl;
	//std::cout << body_angles[0] << std::endl;
	//std::cout << body_angles[1] << std::endl;
	//std::cout << body_angles[2] << std::endl;

	Eigen::VectorXd desPosBase(6);
	desPosBase << 0.0, 0.0, desPosSwingFoot(2)/2, OptComPos(0)-0*0.05*body_angles[0], OptComPos(1)-0*0.01*body_angles[1], OptComPos(2)-0*0.05*0.05*body_angles[2];    //solver->getOptimalCoMPosition();

	//std::cout << "traiettoria swing foot" << std::endl;
	//std::cout << desPosSwingFoot << std::endl;

	//std::cout << "traiettoria base" << std::endl;
	//std::cout << desPosBase << std::endl;

    //    std::cout << "getOptimalCoMPosition" << std::endl;
	//std::cout << solver->getOptimalCoMPosition() << std::endl;   // Think of adding the angular momentum controller here...


    //    std::cout << "Com Angular Velocity" << std::endl;
	//std::cout << mBase->getAngularVelocity() << std::endl;   // See what happens here...

    //    std::cout << "Angular Momentum" << std::endl;
	//std::cout << mBase->getAngularMomentum() << std::endl;   // See what happens here...


	// Compute inverse kinematics
	return getJointVelocitiesQp(desPosBase, actPosBase, desPosSwingFoot, actPosSwingFoot);
        //return getJointVelocitiesStacked(desPosBase, actPosBase, desPosSwingFoot, actPosSwingFoot);
}

Eigen::VectorXd Controller::generateBalance(){

	// Retrieve positions and orientations relative to support foot
	Eigen::VectorXd comCurrentPosition;
	if (balancePoint == TORSO) {
		comCurrentPosition = mBase->getCOM(mSupportFoot);
	} else {
		comCurrentPosition = mRobot->getCOM(mSupportFoot);
	}

	Eigen::VectorXd actPosSwingFoot(6);
	Eigen::VectorXd actPosBase(6);
	actPosSwingFoot << getRPY(mSwingFoot, mSupportFoot), mSwingFoot->getCOM(mSupportFoot);
	actPosBase << getRPY(mBase, mSupportFoot), comCurrentPosition;

	Eigen::VectorXd desPosBase(6);
	desPosBase = balanceBasePos;

	Eigen::VectorXd desPosSwingFoot(6);
	desPosSwingFoot = balanceFootPos;
	// std::cout << "traiettoria swing foot" << std::endl;
	// std::cout << desPosSwingFoot << std::endl;

	// std::cout << "traiettoria base" << std::endl;
	// std::cout << desPosBase << std::endl;
	// Compute inverse kinematics

	// Emilian questo si!	
	if(beheavior == BALANCE) {
		Eigen::VectorXd jointVel = getJointVelocitiesQp(desPosBase, actPosBase, desPosSwingFoot, actPosSwingFoot);
		return jointVel;
	} 
	Eigen::VectorXd zeroVel(24);
	zeroVel.setZero();
	return zeroVel;

	// Emilian questo no!
	// return getJointVelocitiesStacked(desPosBase, actPosBase, desPosSwingFoot, actPosSwingFoot);
}

Eigen::MatrixXd Controller::getTorsoAndSwfJacobian(){

	Eigen::MatrixXd Jacobian_supportToBase;

	if (balancePoint == TORSO) {
		Jacobian_supportToBase =  mRobot->getJacobian(mBase,mSupportFoot) - mRobot->getJacobian(mSupportFoot,mSupportFoot);
	} else {
		Jacobian_supportToBase =  (mRobot->getCOMJacobian(mSupportFoot) - mRobot->getJacobian(mSupportFoot,mSupportFoot));

		Jacobian_supportToBase.col(6).setZero();
		Jacobian_supportToBase.col(7).setZero();
		Jacobian_supportToBase.col(14).setZero();
		Jacobian_supportToBase.col(15).setZero();
		Jacobian_supportToBase.col(16).setZero();
		Jacobian_supportToBase.col(17).setZero();
		Jacobian_supportToBase.col(18).setZero();
		Jacobian_supportToBase.col(25).setZero();
		Jacobian_supportToBase.col(26).setZero();
		Jacobian_supportToBase.col(27).setZero();
		Jacobian_supportToBase.col(28).setZero();
		Jacobian_supportToBase.col(29).setZero();
	}

	Eigen::MatrixXd Jacobian_SupportToSwing =  mRobot->getJacobian(mSwingFoot,mSupportFoot) - mRobot->getJacobian(mSupportFoot,mSupportFoot);

	Eigen::MatrixXd Jacobian_tot_(12, 30);

	Jacobian_tot_ << Jacobian_supportToBase, Jacobian_SupportToSwing;

	Eigen::MatrixXd Jacobian_tot(12, 24);

	// Remove the floating base columns
	Jacobian_tot = Jacobian_tot_.block<12,24>(0, 6);

	//std::cout << sqrt((Jacobian_tot*Jacobian_tot.transpose()).determinant()) << std::endl;

	return Jacobian_tot;
}

Eigen::VectorXd Controller::getJointVelocitiesQp(Eigen::VectorXd desider_pos_base, Eigen::VectorXd actPosBase,
		Eigen::VectorXd desider_pos_SwingFoot, Eigen::VectorXd actPosSwingFoot){

	int nVariables = 24;

	Eigen::VectorXd refPosture = initialConfiguration.segment(6,24);
	Eigen::VectorXd currentPosture = mRobot->getPositions().segment(6,24);
	//std::cout << "redPosture" << refPosture << std::endl;
	//std::cout << "currentPose" << currentPosture << std::endl;
	double jointVelocitiesGain = 0.00000001;
	double postureGain = 0;
	Eigen::MatrixXd taskGain = Eigen::MatrixXd::Identity(12,12);

	// Torso Orientation
	taskGain(0,0) = 1;
	taskGain(1,1) = 1;//0.001;
	taskGain(2,2) = 0;//0.001;

	// CoM Position
	taskGain(5,5) = 1;

	// Swing Foot Orientation
	taskGain(6,6) = 1;
	taskGain(7,7) = 1;
	taskGain(8,8) = 1;

	// Swing Foot Position
	taskGain(11,11) = 1;

	//std::cout << taskGain(4,4) << std::endl;

	Eigen::VectorXd desired_pos(12);
	desired_pos << desider_pos_base, desider_pos_SwingFoot;
	Eigen::VectorXd actual_pos(12);
	actual_pos << actPosBase, actPosSwingFoot;
    //std::cout << "voglio" << desired_pos << std::endl;
	//std::cout << "sono" << actual_pos << std::endl;
	// Cost Function
	Eigen::MatrixXd Jacobian_tot = getTorsoAndSwfJacobian();
	Eigen::MatrixXd costFunctionH = mWorld->getTimeStep()*mWorld->getTimeStep()*Jacobian_tot.transpose()*taskGain*Jacobian_tot +
			postureGain*mWorld->getTimeStep()*mWorld->getTimeStep()*Eigen::MatrixXd::Identity(nVariables,nVariables) +
			jointVelocitiesGain*Eigen::MatrixXd::Identity(nVariables,nVariables);
	//std::cout << "H" << costFunctionH << std::endl;
	Eigen::VectorXd costFunctionF = mWorld->getTimeStep()*Jacobian_tot.transpose()*taskGain*(actual_pos - desired_pos) + postureGain*mWorld->getTimeStep()*(currentPosture - refPosture);
	//std::cout << "F" << costFunctionF << std::endl;
	// Constraint RHipYawPitch and LHipYawPitch to be at the same angle
	Eigen::MatrixXd AHip = Eigen::MatrixXd::Zero(1,nVariables);
	AHip(0,2) = -mWorld->getTimeStep();
	AHip(0,13) = mWorld->getTimeStep();
	Eigen::VectorXd bHip(1);
	bHip(0) = mRobot->getPosition(8) - mRobot->getPosition(19);

	// Constraint max joint accelerations
	double maxAcc = 2;
	Eigen::MatrixXd AAcc = Eigen::MatrixXd::Identity(nVariables,nVariables);
	Eigen::VectorXd bAccMax = 0*qDotOld + Eigen::VectorXd::Ones(nVariables)*maxAcc;
	Eigen::VectorXd bAccMin = 0*qDotOld - Eigen::VectorXd::Ones(nVariables)*maxAcc;

	// Stack Constraints
	//Eigen::MatrixXd A;
	//Eigen::VectorXd bMin;
	//Eigen::VectorXd bMax;
	//A.resize(1+nVariables,nVariables);
	//bMin.resize(1+nVariables);
	//bMax.resize(1+nVariables);
	//std::cout << "bo" << std::endl;
	//A << AHip, AAcc;
	//bMin << bHip, bAccMin;
	//bMax << bHip, bAccMax;

	// Solve the QP
	Eigen::VectorXd solution = solveQP(costFunctionH, costFunctionF, AHip, bHip, bHip);
	return solution;
}

Eigen::VectorXd Controller::getJointVelocitiesStacked(Eigen::VectorXd desider_pos_base, Eigen::VectorXd actPosBase,
		Eigen::VectorXd desider_pos_SwingFoot, Eigen::VectorXd actPosSwingFoot){

	Eigen::VectorXd desired_pos(12);
	desired_pos << desider_pos_base, desider_pos_SwingFoot;

	// Assemble actual positions and orientations
	Eigen::VectorXd actual_pos(12);
	actual_pos << actPosBase, actPosSwingFoot;

	// Get the proper jacobian and pseudoinvert it
	Eigen::MatrixXd Jacobian_tot = getTorsoAndSwfJacobian();
	Eigen::MatrixXd PseudoJacobian_tot = (Jacobian_tot.transpose())*(Jacobian_tot*Jacobian_tot.transpose()).inverse();

	Eigen::VectorXd qDot(24);
	qDot = PseudoJacobian_tot*(ikGain*(desired_pos - actual_pos));

	return qDot;
}

Eigen::VectorXd Controller::getJointVelocitiesStacked(Eigen::VectorXd desider_vel_base, Eigen::VectorXd actVelBase,
		Eigen::VectorXd desider_pos_base, Eigen::VectorXd actPosBase,
		Eigen::VectorXd desider_vel_SwingFoot, Eigen::VectorXd actVelSwingFoot,
		Eigen::VectorXd desider_pos_SwingFoot, Eigen::VectorXd actPosSwingFoot){

	Eigen::VectorXd desired_pos(12);
	Eigen::VectorXd desired_vel(12);
	desired_pos << desider_pos_base, desider_pos_SwingFoot;
	desired_vel << desider_vel_base, desider_vel_SwingFoot;

	// Assemble actual positions and orientations
	Eigen::VectorXd actual_pos(12);
	actual_pos << actPosBase, actPosSwingFoot;

	// Get the proper jacobian and pseudoinvert it
	Eigen::MatrixXd Jacobian_tot = getTorsoAndSwfJacobian();
	Eigen::MatrixXd PseudoJacobian_tot = (Jacobian_tot.transpose())*(Jacobian_tot*Jacobian_tot.transpose()).inverse();

	Eigen::VectorXd qDot(24);
	qDot = PseudoJacobian_tot*(desired_vel + ikGain*(desired_pos - actual_pos));

	return qDot;
}

Eigen::VectorXd Controller::getJointVelocitiesPrioritized(Eigen::VectorXd desider_vel_base, Eigen::VectorXd actVelBase,
		Eigen::VectorXd desider_pos_base, Eigen::VectorXd actPosBase,
		Eigen::VectorXd desider_vel_SwingFoot, Eigen::VectorXd actVelSwingFoot,
		Eigen::VectorXd desider_pos_SwingFoot, Eigen::VectorXd actPosSwingFoot){

	Eigen::MatrixXd JacobianStS_ = -(mRobot->getJacobian(mSupportFoot) - mRobot->getJacobian(mSwingFoot));
	Eigen::MatrixXd JacobianStB_ = (mRobot->getJacobian(mRobot->getBodyNode("base_link")) - mRobot->getJacobian(mSupportFoot));
	Eigen::MatrixXd JacobianStS = JacobianStS_.block<6,24>(0,6);
	Eigen::MatrixXd JacobianStB = JacobianStB_.block<6,24>(0,6);

	Eigen::MatrixXd PinvStS = (JacobianStS.transpose())*(JacobianStS*JacobianStS.transpose()).inverse();
	Eigen::MatrixXd PinvStB = (JacobianStB.transpose())*(JacobianStB*JacobianStB.transpose()).inverse();

	Eigen::VectorXd qDot_base(24);
	Eigen::VectorXd qDot(24);

	double KP = 10;
	double KD = 0.5;

	qDot_base = PinvStB*(desider_vel_base+KP*(desider_pos_base-actPosBase)+KD*(desider_vel_base-actVelBase));
	qDot = PinvStS*(desider_vel_SwingFoot+KP*(desider_pos_SwingFoot-actPosSwingFoot)+KD*(desider_vel_SwingFoot-actVelSwingFoot))
					+(Eigen::MatrixXd::Identity(24,24)-PinvStS*JacobianStS)*qDot_base;

	return qDot;
}

Eigen::Vector3d Controller::getRPY(dart::dynamics::BodyNode* body, dart::dynamics::BodyNode* referenceFrame) {

	Eigen::MatrixXd rotMatrix = body->getTransform(referenceFrame).rotation();

	Eigen::Vector3d RPY;
	RPY << atan2(rotMatrix(2,1),rotMatrix(2,2)),
			atan2(-rotMatrix(2,0),sqrt(rotMatrix(2,1)*rotMatrix(2,1)+rotMatrix(2,2)*rotMatrix(2,2))),
			atan2(rotMatrix(1,0),rotMatrix(0,0));

	return RPY;
	//return body->getTransform(referenceFrame).rotation().eulerAngles(0,1,2);
}

std::pair<Eigen::Vector3d, Eigen::Vector2d> Controller::getZmpFromExternalForces()
{
	// external forces block

	std::vector<double> zmp_v;
	bool left_contact=false;
	bool right_contact =false;
	dart::dynamics::BodyNode* left_foot = mRobot->getBodyNode("l_sole");
	Eigen::Vector3d left_cop;

	// std::cout << "impluso L" << left_foot->getConstraintImpulse() << std::endl;
	if(abs(left_foot->getConstraintImpulse()[5])>0.01){
		left_cop << -(left_foot)->getConstraintImpulse()(1)/(left_foot)->getConstraintImpulse()(5),  (left_foot)->getConstraintImpulse()(0)/(left_foot)->getConstraintImpulse()(5),0.0;
		Eigen::Matrix3d iRotation = left_foot->getWorldTransform().rotation();
		Eigen::Vector3d iTransl   = left_foot->getWorldTransform().translation();
		left_cop = iTransl + iRotation*left_cop;
		left_contact = true;
	}
	dart::dynamics::BodyNode* right_foot = mRobot->getBodyNode("r_sole");
	Eigen::Vector3d right_cop;

	// std::cout << "impluso R" << right_foot->getConstraintImpulse() << std::endl;
	if(abs(right_foot->getConstraintImpulse()[5])>0.01){
		right_cop << -(right_foot)->getConstraintImpulse()(1)/(right_foot)->getConstraintImpulse()(5),(right_foot)->getConstraintImpulse()(0)/(right_foot)->getConstraintImpulse()(5),0.0;
		Eigen::Matrix3d iRotation = right_foot->getWorldTransform().rotation();
		Eigen::Vector3d iTransl   = right_foot->getWorldTransform().translation();
		right_cop = iTransl + iRotation*right_cop;
		right_contact = true;
	}

	if(left_contact && right_contact){
		zmp_v.push_back((left_cop(0)*left_foot->getConstraintImpulse()[5] + right_cop(0)*right_foot->getConstraintImpulse()[5])/(left_foot->getConstraintImpulse()[5] + right_foot->getConstraintImpulse()[5]) );
		zmp_v.push_back((left_cop(1)*left_foot->getConstraintImpulse()[5] + right_cop(1)*right_foot->getConstraintImpulse()[5])/(left_foot->getConstraintImpulse()[5] + right_foot->getConstraintImpulse()[5]) );
		zmp_v.push_back(0.0);
	}else if(left_contact){
		zmp_v.push_back(left_cop(0));
		zmp_v.push_back(left_cop(1));
		zmp_v.push_back(0.0);
	}else if(right_contact){
		zmp_v.push_back(right_cop(0));
		zmp_v.push_back(right_cop(1));
		zmp_v.push_back(0.0);
	}else{
		// No contact detected
		zmp_v.push_back(0.0);
		zmp_v.push_back(0.0);
		zmp_v.push_back(0.0);
	}

	Eigen::Vector3d returnZmp;
	Eigen::Vector2d grf;
	grf << left_foot->getConstraintImpulse()[5], right_foot->getConstraintImpulse()[5];
	returnZmp <<  zmp_v[0], zmp_v[1], zmp_v[2];
	return std::make_pair(returnZmp, grf);
}

std::pair<Eigen::Vector3d, double> Controller::getZmpFromWrench(){
	std::vector<double> zmp_est;
	const dart::collision::CollisionResult& col_res = mWorld->getLastCollisionResult();
	std::vector<dart::collision::Contact> contacts = col_res.getContacts();

	int contact_counter = 1;
	double x_zmp =0.0;
	double y_zmp =0.0;
	double z_zmp =0.0;
	double total_fz = 0;
	if (contacts.size() > 0) {
		for (auto cont : contacts) {

			contact_counter = contact_counter + 1;

			double contactForce;
			if (cont.force.z() < 0) contactForce = 0;
			else 
				contactForce = cont.force.z();

			x_zmp = x_zmp + cont.point.x()*contactForce;
			y_zmp = y_zmp + cont.point.y()*contactForce;
			z_zmp = z_zmp + cont.point.z()*contactForce;
			total_fz = total_fz + contactForce;

		}
	}
	if(total_fz!=0){
		x_zmp = x_zmp/total_fz;
		y_zmp = y_zmp/total_fz;
		z_zmp = z_zmp/total_fz;
	}else{
		x_zmp = 0.0;
		y_zmp = 0.0;
		z_zmp = 0.0;
	}

	zmp_est.push_back(x_zmp);
	zmp_est.push_back(y_zmp);
	zmp_est.push_back(z_zmp);

	Eigen::Vector3d streamable;
	streamable << zmp_est[0], zmp_est[1], zmp_est[2];

	return std::make_pair(streamable, total_fz);
}

Eigen::Vector3d Controller::getZmpFromAngularMomentum(){
	Eigen::Vector3d totalMom;
	Eigen::Vector3d realZMP;

	totalMom.resize(3);
	totalMom.setZero();

	realZMP.resize(3);
	realZMP.setZero();

	Eigen::Vector3d CoMPosition = mRobot->getCOM();
	Eigen::Vector3d CoMAcceleration = mRobot->getCOMLinearAcceleration();

	Eigen::VectorXd oldMom = totalMom;

	totalMom.setZero();

	for (unsigned int i=0; i<mRobot->getNumBodyNodes(); i++){
		Eigen::Vector3d LinearMom = (mRobot->getBodyNode(i))->getLinearMomentum();
		Eigen::Vector3d AngularMom = (mRobot->getBodyNode(i))->getAngularMomentum();

		Eigen::Vector3d iCoMPosition = (mRobot->getBodyNode(i))->getCOM();
		Eigen::Matrix3d iRotation = (mRobot->getBodyNode(i))->getWorldTransform().rotation();

		// transform the momentums in world frame
		LinearMom = iRotation*LinearMom;
		AngularMom = iRotation*AngularMom;

		Eigen::Vector3d d2;
		d2 << (iCoMPosition(0) - CoMPosition(0))*(iCoMPosition(0) - CoMPosition(0)),
				(iCoMPosition(1) - CoMPosition(1))*(iCoMPosition(1) - CoMPosition(1)),
				(iCoMPosition(2) - CoMPosition(2))*(iCoMPosition(2) - CoMPosition(2));

		totalMom = totalMom + (iCoMPosition-CoMPosition).cross(LinearMom);// + AngularMom + (mRobot->getBodyNode(i)->getMass())*d2;
	}

	Eigen::VectorXd momDot = (totalMom - oldMom)/0.01;

	realZMP << CoMPosition(0) - (CoMPosition(2)/(9.8-0*CoMAcceleration(2))) * CoMAcceleration(0)
			    		 + (1/(5*(9.8-0*CoMAcceleration(2)))) * momDot(1),
						 CoMPosition(1) - (CoMPosition(2)/(9.8-0*CoMAcceleration(2))) * CoMAcceleration(1)
						 + (1/(5*(9.8-0*CoMAcceleration(2)))) * momDot(0),
						 0;

	return realZMP;
}

dart::dynamics::SkeletonPtr Controller::getRobot() const {
	return mRobot;
}

dart::dynamics::BodyNode* Controller::getSupportFoot() {
	return mSupportFoot;
}

mpcSolver::MPCSolver* Controller::getSolver() {
	return solver;
}

Eigen::VectorXd Controller::getEndEffector() const {
	return mRobot->getCOM();
}

void Controller::keyboard(unsigned char /*_key*/, int /*_x*/, int /*_y*/) {

}


Eigen::VectorXd Controller::getOmnidirectionalSwingFootTrajectoryMPC(Eigen::VectorXd targetFootstepPosition, Eigen::VectorXd swingFootStartingPosition,
		Eigen::VectorXd swingFootActualPosition, double stepHeight, int indInitial){

	double delta = mWorld->getTimeStep();
	int ind = mWorld->getSimFrames();

	double time=ind*delta;
	double ti=indInitial*delta;

	int S = (int)(singleSupportDuration/delta);
	int D = (int)(doubleSupportDuration/delta);

	// Generate XY polynomial for single support
	Eigen::VectorXd boundaryConditions(6);
	Eigen::MatrixXd polyMatrix(6,6);
	double dur = delta*S;

	polyMatrix << 0,0,0,0,0,1,
			pow(dur,5),pow(dur,4),pow(dur,3),pow(dur,2),pow(dur,1),1,
			0,0,0,0,1,0,
			5*pow(dur,4),4*pow(dur,3),3*pow(dur,2),2*pow(dur,1),1,0,
			0,0,0,2,0,0,
			5*4*pow(dur,3),4*3*pow(dur,2),3*2*pow(dur,1),2,0,0;

	boundaryConditions << 0,1,0,0,0,0;
	Eigen::VectorXd coeffs = polyMatrix.inverse()*boundaryConditions;

	double t = (ind-indInitial)*delta;
	double polyInT = coeffs(0)*pow(t,5) + coeffs(1)*pow(t,4) + coeffs(2)*pow(t,3) + coeffs(3)*pow(t,2) + coeffs(4)*pow(t,1) + coeffs(5);
//	double polyDotInT = 5*coeffs(0)*pow(t,4) + 4*coeffs(1)*pow(t,3) + 3*coeffs(2)*pow(t,2) + 2*coeffs(3)*pow(t,1) + coeffs(4);

	Eigen::VectorXd swingFootPosition = Eigen::VectorXd::Zero(6);
	// std::cout << "ind" << ind << std::endl;
	// std::cout << "indInitial" << indInitial << std::endl;
	// std::cout << "t" << t << std::endl;
	// std::cout << "polyMatrix" << polyMatrix << std::endl;
	// std::cout << "polyInt" << polyInT << std::endl;
	
	if (ind<indInitial){
		swingFootPosition(3) = swingFootStartingPosition(3);
		swingFootPosition(4) = swingFootStartingPosition(4);

	} else if (ind>=indInitial+S || (ind>=indInitial+S/2 && swingFootActualPosition(5)<=0) ){
		swingFootPosition(3) = targetFootstepPosition(0);
		swingFootPosition(4) = targetFootstepPosition(1);

		swingFootPosition(2) = targetFootstepPosition(3);

	} else {
		swingFootPosition(3) = swingFootStartingPosition(3) + (targetFootstepPosition(0) - swingFootStartingPosition(3))*polyInT;
		swingFootPosition(4) = swingFootStartingPosition(4) + (targetFootstepPosition(1) - swingFootStartingPosition(4))*polyInT;
		swingFootPosition(5) = -(4*stepHeight/pow((S*delta),2))*(time-ti)*(time-ti-S*delta);

		swingFootPosition(2) = targetFootstepPosition(3)*polyInT;

		/*swingFootPosition(3) = swingFootActualPosition(3) + (targetFootstepPosition(0) - swingFootActualPosition(3))*0.5;
		swingFootPosition(4) = swingFootActualPosition(4) + (targetFootstepPosition(1) - swingFootActualPosition(4))*0.5;
		swingFootPosition(5) = -(4*stepHeight/pow((S*delta),2))*(time-ti)*(time-ti-S*delta);

		swingFootPosition(2) = targetFootstepPosition(3)*polyInT;*/
	}
	// std::cout << "swingFootPosition" << swingFootPosition << std::endl;

	return swingFootPosition;
}

void Controller::setInitialConfiguration() {

	Eigen::VectorXd floatingBaseConfiguration(6);
	Eigen::VectorXd headConfiguration(2);
	Eigen::VectorXd rightArmConfiguration(5);
	Eigen::VectorXd leftArmConfiguration(5);
	Eigen::VectorXd rightLegConfiguration(6);
	Eigen::VectorXd leftLegConfiguration(6);

	floatingBaseConfiguration << 0.0, 3*M_PI/180, 0.0, 0.0, 0.0, 0.35;
	headConfiguration << 0.0, 0.0;
	rightArmConfiguration << 80*M_PI/180, -10*M_PI/180.0, 50*M_PI/180.0, 2*M_PI/180.0, 0.0;
	leftArmConfiguration << 80*M_PI/180, 10*M_PI/180.0, -50*M_PI/180.0, -2*M_PI/180.0, 0.0;
	rightLegConfiguration << 0.0, 0.0, -20*M_PI/180.0, 34*M_PI/180.0, -19*M_PI/180.0, 0.0;
	leftLegConfiguration << 0.0, 0.0, -20*M_PI/180.0, 34*M_PI/180.0, -19*M_PI/180.0, 0.0;

	initialConfiguration << floatingBaseConfiguration, headConfiguration, leftLegConfiguration, leftArmConfiguration, rightLegConfiguration, rightArmConfiguration;

	mRobot->setPositions(initialConfiguration);
}

void Controller::setComTargetHeight(double height) {
	solver->setComTargetHeight(height);
	// TODO: change target height in observer
}

void Controller::setReferenceVelocityX(double ref) {
	solver->setReferenceVelocityX(ref);
}

void Controller::setReferenceVelocityY(double ref) {
	solver->setReferenceVelocityY(ref);
}

void Controller::setReferenceVelocityOmega(double ref) {
	solver->setReferenceVelocityOmega(ref);
}

void Controller::setBalancePointCom() {
	balancePoint = COM;
}

void Controller::setBalancePointTorso() {
	balancePoint = TORSO;
}

void Controller::setBeheaviorBalance() {
	beheavior = BALANCE;
}

void Controller::setBeheaviorWalk() {
	beheavior = WALK;
}

void Controller::setBeheaviorFreeBalance() {
	beheavior = FREE_BALANCE;
}

void Controller::setVelGain(double value) {
	solver->setVelGain(value);
}

void Controller::setZmpGain(double value) {
	solver->setZmpGain(value);
}

Eigen::VectorXd Controller::getBalanceBasePos() {
	return balanceBasePos;
}

Eigen::VectorXd Controller::getBalanceFootPos() {
	return balanceFootPos;
}

void Controller::setBalanceBasePos(Eigen::VectorXd pos) {
	balanceBasePos = pos;
}

void Controller::setBalanceFootPos(Eigen::VectorXd pos) {
	balanceFootPos = pos;
}


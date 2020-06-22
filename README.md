-- For Emilian

Install https://launchpad.net/ubuntu/+archive/primary/+sourcefiles/openscenegraph/3.2.3+dfsg1-3ubuntu2/openscenegraph_3.2.3+dfsg1.orig.tar.gz and
and https://github.com/assimp/assimp/releases/tag/v4.1.0/ from source. Then add paths to LD_LIBRARY_PATH. Then use the modified version of CMakeLists.txt and afterwards before make, comment the bad imports.


# Autonomous and Mobile Robotics 2020
Project repository for the course of Autonomous and Mobile Robotics, Sapienza University of Rome.

##  External force observers for the NAO
### Deployment of a wrapper for external force disturbance observers acting on a NAO robot while walking

This repository contains the C++ implementation of a Luenberger, a Kalman-Filter based and a Stephens' observers for estimating an unknown external disturbance acting on the NAO's center of mass. Moreover a C++ wrapper for an easier interface with the observers and a script to plot and save the results of the observers were implemented.


Contents:

- [Structure](#structure)
- [Execution](#execution)
- [Plot Results](#plot-results)
- [References](#references)
- [Team Members](#team-members)

## Structure
- To be completed.

## Execution
1. Install DART following https://dartsim.github.io/install_dart_on_ubuntu.html
2. Go into intrinsically_stable_mpc and create a build directory:
```
mkdir build
cd build
```
3. Build and compile the simulation:
```
cmake ..
make
```
4. Use in-app widget to customize the experiments.

## Plot Results
Using `plot.py` it is possible to plot and save the logs produced by `controller.cpp`. To run the script specify the following parameters:
 * `--pathESTx`: (optional) path to the .txt file of the estimates of the observer along the x axis
 * `--pathGTx`: (optional) path to the .txt file of the ground truth values along the x axis
 * `--pathESTy`: (optional) path to the .txt file of the estimates of the observer along the y axis
 * `--pathGTy`: (optional) path to the .txt file of the ground truth values along the y axis
 * `--pathESTz`: (optional) path to the .txt file of the estimates of the observer along the z axis
 * `--pathGTz`: (optional) path to the .txt file of the ground truth values along the z axis
 * `--pathSavings`: path to the observer folder where to save the plots (default `./plots/`)
 * `--margin`: float value for the confidence interval of the ground truth (default `1.0`)
 * `--observer`: name of the observer (default `./plots/`)
 * `--quantity`: (optional) quantity to plot, if `None`, all quantities will be plotted.
 
 
 More info available by calling `--h`. An example is 

```
python plot.py --pathESTx ./logs/stephens_x_est.txt --pathGTx ./logs/stephens_x_gt.txt --observer stephens --quantity force
```

## References

The papers used in this work include, but are not limited to:
	
- N. Scianca, M. Cognetti, D. De Simone, L. Lanari, G. Oriolo. *Intrinsically stable MPC for humanoid gait generation.* 	2016 IEEE-RAS 16th International Conference on Humanoid Robots (Humanoids), Cancun, 2016, pp. 601-606, doi: 10.1109/HUMANOIDS.2016.7803336.
- F. M. Smaldone, N. Scianca, V. Modugno, L. Lanari, G. Oriolo. *Gait Generation using Intrinsically Stable MPC in the Presence of Persistent Disturbances.* 2019 IEEE-RAS 19th International Conference on Humanoid Robots (Humanoids), Toronto, ON, Canada, 2019, pp. 651-656, doi: 10.1109/Humanoids43949.2019.9035068.
- L. Hawley, W. Suleiman. *External force observer for medium-sized humanoid robots.*	2016 IEEE-RAS 16th International Conference on Humanoid Robots (Humanoids), Cancun, 2016, pp. 366-371, doi: 10.1109/HUMANOIDS.2016.7803302.
- L. Hawley, R. Rahem, W. Suleiman. *Kalman Filter Based Observer for an External Force Applied to Medium-sized Humanoid Robots.* 2018 IEEE/RSJ International Conference on Intelligent Robots and Systems (IROS), Madrid, 2018, pp. 1204-1211, doi: 10.1109/IROS.2018.8593610.
- B. J. Stephens. *State estimation for force-controlled humanoid balance using simple models in the presence of modeling error.* 2011 IEEE International Conference on Robotics and Automation, Shanghai, 2011, pp. 3994-3999, doi: 10.1109/ICRA.2011.5980358.

## Team members

- Postolache Emilian, 1649271.
- Ratini Riccardo, 1656801.
- Simionato Giada, 1822614.

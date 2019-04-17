# Sensor Fusion Assignment - Autonomous Vehicle Dynamics and Control


### Discription
This Software Project was created as part of the Sensor Fusion module which was held by Dr. Hyo-Sang Shin, Cranfield University for the students of the master course Autonomous Vehicle Dynamics and Control. This repository contains MATLAB software for the Performance Analsysis of multiple Sensor Fusion Algorithms for the Position Estimation of a non-moving GPS Jammer from a Unmanned Aerial Vehicle Platform. 

The obstacle of the assignment was the development and implementation of different sensor fusion algorithms in a synthetic environment provided in MATLAB/Simulink. The provided problem definition was divided into two sub-parts. The first part is focusing on the development of certain sensor fusion algorithms in an environment with isotropic jammer behaviour. This section is concentrating on the implementation of the Extended Kalman, Unscented Kalman and Particle Filter algorithms to estimate the position of a none moving GPS jamming platform. In the second part, the jammer patterns are anisotropic in theire spreading
behaviour. The prior developed algorithms (Extended Kalman, Unscented Kalman and Particle Filters) should be used to analyse their behaviour with the new jammer patterns. As last part of the assignment, new abd unpresented Sensor Fusion algorithms should have been implemented. The goal for that algorithms was to improve the position estimation for the case of having an anisotropic jammer behaviour

#### GPS Jammer Localisation Methods

The position of the GPS jamming vehicle is not directly observable and must be therefore observed by utilizing obseravble measurements. In order to localise the target that is confusing the GPS signals, a power measurement of the jamming signal is used. This method is called the Received signal strength (RSS) method, which is often used for localisation of energy emitting sources.

#### UAV Guidance Approach

The guidance approach which is utilised for the UAV guidance is a vector field based path following. In which the UAV first follow a straight line approach to the estimated position of the target. When the UAV reaches a certain distance it starts to follow a loiter path.


#### Sensor Fusion Algorithms

The goal of the assignment was it to develop and implement different sensor fusion algorihtms for both an isotropic and anistropic GPS jammer pattern. The following algorhtms were implemented:
* Extended Kalman Filter
* Unscendet Kalman Filter
* Particle Filter
* Extended Particle Filter
* Unscendet Particle Filter
* H-Infinity Filter
* Adaptive Kalman Filter
* H-Infinity Particle Filter


#### Particle Filter Resampling

A set of potential resampling approaches have been implmented. The goal of that analysis was to compare the different performances of the resampling methods in the same environment and to choose the algorithm with the best performance. Following methods have been implemented and investigated:

* Multinomial Resampling
* Systemematic Resampling
* Residual Resampling
* Residual Systematic Resampling
* Local Selection Resampling
* Stratified Resampling 

### Usage

the isotropic and the anisotropic algorithms are both in seperated directories. In order to get an visualised test you can run the Main_isotropic_XXX or Main_anisotropic_XXX files. XXX stands for the name of the sensor fusion algorithm which is used for that test scenario.


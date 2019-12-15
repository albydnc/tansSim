# TANS 2019-2020 Final Project
### Description
This program is a basic implementation of simulation and recontruction of particles colliding events.


### Simulation
The program simulates a detector made of a beampipe and two cylindrical layers.
Collision vertexes and multiplicity are generated on a normal distribution and particles follow the pseudorapidity and angle distribution from a given kinematic file.

##### Interaction with detector layers
The interaction with detectors is following a basic linear intesection.
Gaussian smearing is applied to the points and multiple scattering can be enabled.
It is also possible to generate a fixed number of random noise points.
Data is saved in `simulation.root` file.

### Reconstruction
The recontruction part reads the events from the previous simulation file to reconstruct events' vertexes.

##### Tracklet Algorithm
A simple linear interpolation is done after a large angle filter. The results are then plotted and the mode value is returned.

##### Statistics generation
The reconstruction process will print the resulting statistics, giving the user histograms on:
+ Vertex resolution distribution
+ Resolution vs Z-axis span
+ Efficiency vs Multiplicity
+ Efficiency vs Z-axis span

### Simulation and Reconstruction Performance
**Testbench:**
+ Ubuntu 18.04 LTS
+ ROOT v6.18
+ Intel i7-6700HQ Laptop CPU
+ 16GB DDR4 RAM
+ 512GB NVMe SSD

Simulation of 1M events takes 80 seconds on average.
Reconstruction of 1M events takes 40 seconds on average.

Performance is greatly improved by using the C++ STL data classes over _TClonesArray_.
Simulation of 1M events is around 10 seconds average. **8x faster**.
As for now, _TTree_ does not support `std::vector<std::array<double,2>>` so the development is not further going on.
I will keep track of ROOT new releases to see if it will be supported.

*Written by Alberto Perro - 2019*

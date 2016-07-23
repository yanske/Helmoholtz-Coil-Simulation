#include <iostream>
#include <fstream>
#include "coil.h"
#include "bfield.h"

/*
Simulation of magnetic field caused by two Helmoholtz coils (non-parallel)
Yan Ke
*/

void runSimulation()
{
    std::cout<<"Starting Helmoholtz coil simulation..."<<std::endl;
    double fieldsize = 0.20; //defines size of simulation
    double scale = 0.005; //m per increment
    int dimension = fieldsize/scale; //dimension of array
    Bfield one(dimension, dimension, dimension);
    Bfield two(dimension, dimension, dimension);

    //magnetic field calculated by the integration of the Biot-Savart Law (dB = uidsr/4pi r^3)
    //, where r =P-R, P = center of coil, R = (Rsintheta, 0, -Rcostheta)
    one.calculateBfield(scale, 0.01);
    std::cout<<"Coil 1 simulation done."<<std::endl;
    two.calculateBfield(scale, 0.11);
    std::cout<<"Coil 2 simulation done."<<std::endl;
    //rotate the 2nd field using the rotation matrix
    two.rotateBfield(25); //angle in degrees
    //add the Bfield components
    one.addBfield(two);
    std::cout<<one;

}
int main()
{
    runSimulation();
    return 0;
}

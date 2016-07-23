#ifndef BFIELD_H
#define BFIELD_H
#include "coil.h"
#include <cmath>
#include <math.h>
#include <fstream>
#include <iostream>

const double pi = 3.14159265358979;
class Bfield: public Coil
{
private:
    Point3d m_dimension;
    double *m_data;
    double m_scale;

public:
    //creating a 4-d array with a set 4th dimension
    //x, y, z coordinates for where the magnetic field is, and 3d components of the magnetic field in the position
    Bfield(int x, int y, int z)
    {
        //dimensions of our simulaiton
        m_dimension.x=x;
        m_dimension.y=y;
        m_dimension.z = z;
        m_dimension.length = (x*y*z*3);
        m_data = new double [m_dimension.length];
    }

    //function to return value from a 4d input
    double& getValue(int x, int y, int z, int b)
    {
        return m_data[((x*m_dimension.y+y)*m_dimension.z+z)*3+b];
    }
    void calculateBfield(const double scale, const double initialx)
    {
        m_scale = scale;
        int steps = (m_thickness/m_scale) +0.5;
        double integralstep = 0.1;
        int dn = m_rotation/steps;
        const double mew = 1.2566370614*pow(10,-6);
        for (int x =0; x < m_dimension.x; x++)
        {
            for (int y = 0; y < m_dimension.y; y++)
            {
                for (int z =0; z < m_dimension.z; z++)
                {
                    for (int n = 0; n < steps; n++)
                    {
                        //setup variables here (unique to this coil)
                        //position
                        m_position.x = initialx + n*(m_thickness/steps);
                        m_position.y = m_radius;
                        m_position.z = m_radius;

                        //integrate using Riemann Sum
                        for (double theta = 0.0; theta < 2*pi ; theta += integralstep)
                        {
                            double expressionconstant = (dn)*(integralstep)*(mew)*(m_io)/(4*pi*(pow((pow((x*scale-m_position.x)-m_radius*sin(theta),2)+pow((y*scale-m_position.y),2)+pow((z*scale-m_position.z)+cos(theta),2)),1.5)));
                            getValue(x,y,z,0) += expressionconstant*(-1*m_radius*sin(theta)*(y*scale-m_position.y));
                            getValue(x,y,z,1) += expressionconstant*(m_radius*sin(theta)*(x*scale-m_position.x)-(pow(m_radius*sin(theta),2))-m_radius*cos(theta)*(z*scale-m_position.z)-(pow(m_radius*cos(theta),2)));
                            getValue(x,y,z,2) += expressionconstant*(m_radius*cos(theta)*(y*scale-m_position.y));
                        }
                    }
                }
            }
        }
    }
    void rotateBfield(const double rotatez) //rotation independent of the z-axis
    {
        double alpha = rotatez*pi/180; //convert angle from radian to degree
        double *dup_data = new double[m_dimension.length];
        double rotationmatrix[4]={ cos(alpha), -sin(alpha), sin(alpha), cos(alpha)};  //rotationmatrix  U= RX
        int u,v; //rotated coordinates
        for (int x =0; x < m_dimension.x; x++)
        {
            for (int y = 0; y < m_dimension.y; y++)
            {
                for (int z =0; z < m_dimension.z; z++)
                {
                    u = rotationmatrix[0]*x+rotationmatrix[1]*y;
                    v = rotationmatrix[2]*x+rotationmatrix[3]*y;
                    //shift data into new pointer, delete old pointer, copy data over
                    if(u>=0 && u<m_dimension.x && v>=0 && v<m_dimension.y)
                    {
                        dup_data[((x*m_dimension.y+y)*m_dimension.z+z)*3+0] = getValue(u,v,z,0);
                        dup_data[((x*m_dimension.y+y)*m_dimension.z+z)*3+1] = getValue(u,v,z,1);
                        dup_data[((x*m_dimension.y+y)*m_dimension.z+z)*3+2] = getValue(u,v,z,2);
                    }
                }
            }
        }
        delete[] m_data;
        m_data = dup_data;
        std::cout<<"Rotation done."<<std::endl;
    }
    //sum values from other magnetic field
    void addBfield(Bfield &otherfield)
    {
        for (int x =0; x < m_dimension.x; x++)
        {
            for (int y = 0; y < m_dimension.y; y++)
            {
                for (int z =0; z < m_dimension.z; z++)
                {
                    getValue(x,y,z,0) += otherfield.getValue(x,y,z,0);
                    getValue(x,y,z,1) += otherfield.getValue(x,y,z,1);
                    getValue(x,y,z,2) += otherfield.getValue(x,y,z,2);
                }
            }
        }
        std::cout<<"Magnetic field combined."<<std::endl;
    }
    //friend function to export Bfield data as CSV file
    friend std::ostream& operator<<(std::ostream &out, Bfield &bfield);
    ~Bfield()
    {
        delete[] m_data;
    }
};

std::ostream& operator<<(std::ostream &out, Bfield &bfield)
{
    std::ofstream fileout("SimulationData.csv");
    if(!fileout)
    {
        std::cerr<<"File error"<<std::endl;
    }
    std::cout<<"Exporting data..."<<std::endl;
    fileout<<"Px,Py,Pz,Bx,By,Bz"<<std::endl;
    for (int x =0; x < bfield.m_dimension.x; x++)
    {
        for (int y = 0; y<bfield.m_dimension.y; y++)
        {
            for (int z =0; z < bfield.m_dimension.z; z++)
            {
               fileout<<x*bfield.m_scale<<","<<y*bfield.m_scale<<","<<z*bfield.m_scale<<","<<bfield.getValue(x,y,z,0)<<","<<bfield.getValue(x,y,z,1)<<","<<bfield.getValue(x,y,z,2)<<std::endl;
            }
        }
    }
    std::cout<<"Simulation Completed."<<std::endl;
}
#endif // BFIELD_H

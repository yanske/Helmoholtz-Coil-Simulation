#ifndef COIL_H
#define COIL_H
struct Point3d
{
    int x;
    int y;
    int z;
    int length;
};

class Coil
{
protected:
    int m_rotation;
    double m_thickness;
    double m_radius;
    double m_io; //current
    Point3d m_position;

public:
    Coil()
        :m_rotation(320), m_thickness(0.0193), m_radius(0.068), m_io(0.340)
    {
        //position of the coil in the defined dimensions (set to 0 for now)
        m_position.x = 0;
        m_position.y = 0;
        m_position.z = 0;
    }
    virtual ~Coil()
    {
    }
};

#endif // COIL_H

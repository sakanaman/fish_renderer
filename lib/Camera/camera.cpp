#include "camera.hpp"
#include <cmath>
#include <iostream>
//utility functions
template<class Real>
void cross(Real* cross, const Real* v, const Real* w)
{
    cross[0] = v[1] * w[2] - v[2] * w[1];
    cross[1] = v[2] * w[0] - v[0] * w[2];
    cross[2] = v[0] * w[1] - v[1] * w[0];
}

template<class Real>
void normalize(const Real* v, Real* result)
{
    Real length = v[0] * v[0] 
                + v[1] * v[1] 
                + v[2] * v[2];

    length = std::sqrt(length);

    result[0] = v[0]/length;
    result[1] = v[1]/length;
    result[2] = v[2]/length;
}

//camera implementation
template<class Real>
void CreateRightVec(const Real* camforward, Real* result){
    Real tmp[3];
    Real base[] = {0.0, 1.0, 0.0};
    cross(tmp, base, camforward);
    normalize(tmp, result);
    
    result[0] = -1 * result[0];
    result[1] = -1 * result[1];
    result[2] = -1 * result[2];
}

template<class Real>
void CreateUpVec(const Real* camForward, const Real* camRight, Real* result){
    Real tmp[3];
    cross(tmp, camForward, camRight);
    normalize(tmp, result);
}

template<class Real>
PinholeCamera<Real>::PinholeCamera(const Real* camP, const Real* camF)
{
    for(int i = 0; i < 3; ++i)
    {
        camPos[i] = camP[i];
        camForward[i] = camF[i];
    }
    CreateRightVec(camForward, camRight);
    CreateUpVec(camForward, camRight, camUp);
}

template<class Real>
void PinholeCamera<Real>::CreateFirstRay(const Real u, const Real v, Real* ray_origin, Real* ray_dir) const
{
    Real sensorPos[3];
    sensorPos[0] = camPos[0] + u * camRight[0] - v * camUp[0];
    sensorPos[1] = camPos[1] + u * camRight[1] - v * camUp[1];
    sensorPos[2] = camPos[2] + u * camRight[2] - v * camUp[2];
    
    Real Pinhole[3];
    Pinhole[0] = camPos[0] + camForward[0];
    Pinhole[1] = camPos[1] + camForward[1];
    Pinhole[2] = camPos[2] + camForward[2];

    ray_origin[0] = sensorPos[0];
    ray_origin[1] = sensorPos[1];
    ray_origin[2] = sensorPos[2];

    // ray_dir[0] = Pinhole[0] - sensorPos[0];
    // ray_dir[1] = Pinhole[1] - sensorPos[1];
    // ray_dir[2] = Pinhole[2] - sensorPos[2];

    Real result[3];
    result[0] = Pinhole[0] - sensorPos[0];
    result[1] = Pinhole[1] - sensorPos[1];
    result[2] = Pinhole[2] - sensorPos[2];
    normalize(result, ray_dir);
}


// template<class Real>
// Real PinholeCamera<Real>::Sensitivity(const Real u, const Real v) const {
//     return u ;
// }

template class PinholeCamera<float>;
template class PinholeCamera<double>;
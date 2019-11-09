#ifndef CAMERA_HPP
#define CAMERA_HPP


//if needed standardization, use template!!
// template<class Real, class CameraDetail>
// class Camera
// {
// public:
//     Camera(CameraDetail cam);
//     void CreateInitRay(Real* origin, Real* dir) const;
// private:
//     CameraDetail camdetail;
// };

template<class Real>
class PinholeCamera
{
public:
    PinholeCamera(const Real* camP, const Real* camF);
    void CreateFirstRay(const Real u, const Real v, Real* ray_origin, Real* ray_dir) const;
    Real Sensitivity(const Real u, const Real v) const;
private:
    Real camPos[3];
    Real camForward[3];
    Real camRight[3];
    Real camUp[3];
};

#endif